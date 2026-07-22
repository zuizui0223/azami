#!/usr/bin/env python3
"""Local blinded bounding-box annotation app for the independent detector audit.

The app uses Tkinter so no browser plug-in or hosted service is required. It loads
only the blinded manifest and downloaded source images. Taxonomy, coordinates and
frozen detector predictions are never displayed.
"""
from __future__ import annotations

import argparse
import csv
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
from PIL import Image, ImageTk

FIELDS = [
    "audit_id", "annotator_id", "annotation_status", "assessability",
    "image_quality", "image_width", "image_height", "gt_index", "gt_label",
    "life_stage", "occlusion", "edge_truncated", "gt_x1", "gt_y1", "gt_x2",
    "gt_y2", "notes",
]
ASSESSABILITY = ("assessable", "no_target", "unassessable")
IMAGE_QUALITY = ("good", "moderate", "poor")
LIFE_STAGE = ("anthesis", "bud", "postanthesis", "uncertain")
OCCLUSION = ("none", "partial", "severe")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packet-dir", required=True)
    parser.add_argument("--annotator-id", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--canvas-width", type=int, default=1100)
    parser.add_argument("--canvas-height", type=int, default=760)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


@dataclass
class BoxRecord:
    x1: float
    y1: float
    x2: float
    y2: float
    life_stage: str
    occlusion: str
    edge_truncated: bool


class AnnotationApp:
    def __init__(self, root: Any, args: argparse.Namespace) -> None:
        import tkinter as tk
        from tkinter import messagebox, ttk

        self.tk = tk
        self.ttk = ttk
        self.messagebox = messagebox
        self.root = root
        self.args = args
        self.packet_dir = Path(args.packet_dir).resolve()
        self.output = Path(args.output).resolve()
        self.annotator_id = args.annotator_id.strip()
        if not self.annotator_id:
            raise ValueError("annotator-id cannot be empty")

        self.manifest = self._find_one("detector_independent_audit_blinded_manifest.csv")
        self.assignments_path = self._find_one("detector_independent_audit_assignments.csv")
        manifest = pd.read_csv(self.manifest, dtype=str, keep_default_na=False)
        assignments = pd.read_csv(self.assignments_path, dtype=str, keep_default_na=False)
        assigned_ids = assignments.loc[
            assignments["annotator_id"].eq(self.annotator_id), "audit_id"
        ].astype(str)
        self.tasks = manifest.loc[manifest["audit_id"].isin(set(assigned_ids))].copy()
        self.tasks = self.tasks.sort_values("audit_id").reset_index(drop=True)
        if self.tasks.empty:
            raise ValueError(f"No tasks assigned to {self.annotator_id}")
        self.saved_rows = self._load_saved()
        self.saved_by_id = self._group_saved(self.saved_rows)
        incomplete = [
            index for index, audit_id in enumerate(self.tasks["audit_id"])
            if audit_id not in self.saved_by_id
        ]
        self.index = incomplete[0] if incomplete else 0
        self.boxes: list[BoxRecord] = []
        self.image: Image.Image | None = None
        self.tk_image: ImageTk.PhotoImage | None = None
        self.scale = 1.0
        self.drag_start: tuple[float, float] | None = None
        self.drag_item: int | None = None

        self.root.title(f"Detector audit annotation — {self.annotator_id}")
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)
        controls = ttk.Frame(root, padding=6)
        controls.pack(side=tk.TOP, fill=tk.X)
        self.progress_label = ttk.Label(controls)
        self.progress_label.grid(row=0, column=0, columnspan=8, sticky="w")

        ttk.Label(controls, text="Status").grid(row=1, column=0)
        self.assessability = tk.StringVar(value="assessable")
        ttk.Combobox(
            controls, textvariable=self.assessability, values=ASSESSABILITY,
            state="readonly", width=14
        ).grid(row=1, column=1)
        ttk.Label(controls, text="Image quality").grid(row=1, column=2)
        self.image_quality = tk.StringVar(value="good")
        ttk.Combobox(
            controls, textvariable=self.image_quality, values=IMAGE_QUALITY,
            state="readonly", width=12
        ).grid(row=1, column=3)
        ttk.Label(controls, text="New box stage").grid(row=1, column=4)
        self.life_stage = tk.StringVar(value="anthesis")
        ttk.Combobox(
            controls, textvariable=self.life_stage, values=LIFE_STAGE,
            state="readonly", width=13
        ).grid(row=1, column=5)
        ttk.Label(controls, text="Occlusion").grid(row=1, column=6)
        self.occlusion = tk.StringVar(value="none")
        ttk.Combobox(
            controls, textvariable=self.occlusion, values=OCCLUSION,
            state="readonly", width=10
        ).grid(row=1, column=7)

        self.edge_truncated = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            controls, text="New box touches image edge", variable=self.edge_truncated
        ).grid(row=2, column=0, columnspan=2, sticky="w")
        ttk.Label(controls, text="Notes").grid(row=2, column=2)
        self.notes = tk.StringVar(value="")
        ttk.Entry(controls, textvariable=self.notes, width=52).grid(
            row=2, column=3, columnspan=5, sticky="ew"
        )

        button_frame = ttk.Frame(root, padding=4)
        button_frame.pack(side=tk.TOP, fill=tk.X)
        ttk.Button(button_frame, text="Previous", command=self.previous).pack(side=tk.LEFT)
        ttk.Button(button_frame, text="Undo box", command=self.undo_box).pack(side=tk.LEFT)
        ttk.Button(button_frame, text="Clear boxes", command=self.clear_boxes).pack(side=tk.LEFT)
        ttk.Button(button_frame, text="Save", command=self.save_current).pack(side=tk.RIGHT)
        ttk.Button(button_frame, text="Save && Next", command=self.save_next).pack(side=tk.RIGHT)

        self.canvas = tk.Canvas(
            root, width=args.canvas_width, height=args.canvas_height,
            background="#202020", cursor="crosshair"
        )
        self.canvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.canvas.bind("<ButtonPress-1>", self.on_press)
        self.canvas.bind("<B1-Motion>", self.on_drag)
        self.canvas.bind("<ButtonRelease-1>", self.on_release)
        self.root.bind("<Control-z>", lambda _: self.undo_box())
        self.root.bind("<Control-s>", lambda _: self.save_current())
        self.root.bind("<Return>", lambda _: self.save_next())
        self.load_task()

    def _find_one(self, filename: str) -> Path:
        matches = list(self.packet_dir.rglob(filename))
        if len(matches) != 1:
            raise FileNotFoundError(f"Expected one {filename}, found {len(matches)}")
        return matches[0]

    def _image_path(self, filename: str) -> Path:
        matches = list(self.packet_dir.rglob(filename))
        if len(matches) != 1:
            raise FileNotFoundError(f"Expected one image {filename}, found {len(matches)}")
        return matches[0]

    def _load_saved(self) -> list[dict[str, Any]]:
        if not self.output.exists() or self.output.stat().st_size == 0:
            return []
        table = pd.read_csv(self.output, dtype=str, keep_default_na=False)
        return table.to_dict("records")

    @staticmethod
    def _group_saved(rows: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
        grouped: dict[str, list[dict[str, Any]]] = {}
        for row in rows:
            grouped.setdefault(text(row.get("audit_id")), []).append(row)
        return grouped

    def load_task(self) -> None:
        row = self.tasks.iloc[self.index]
        audit_id = text(row["audit_id"])
        filename = text(row["screen_download_filename"])
        path = self._image_path(filename)
        self.image = Image.open(path).convert("RGB")
        self.boxes = []
        self.assessability.set("assessable")
        self.image_quality.set("good")
        self.notes.set("")
        saved = self.saved_by_id.get(audit_id, [])
        if saved:
            first = saved[0]
            self.assessability.set(text(first.get("assessability")) or "assessable")
            self.image_quality.set(text(first.get("image_quality")) or "good")
            self.notes.set(text(first.get("notes")))
            for item in saved:
                if text(item.get("gt_index")):
                    self.boxes.append(
                        BoxRecord(
                            float(item["gt_x1"]), float(item["gt_y1"]),
                            float(item["gt_x2"]), float(item["gt_y2"]),
                            text(item.get("life_stage")) or "uncertain",
                            text(item.get("occlusion")) or "none",
                            text(item.get("edge_truncated")).lower() in {"true", "1", "yes"},
                        )
                    )
        self.render()
        complete = len(self.saved_by_id)
        self.progress_label.configure(
            text=f"{self.index + 1}/{len(self.tasks)} — {audit_id} — saved {complete}/{len(self.tasks)}"
        )

    def render(self) -> None:
        if self.image is None:
            return
        from PIL import ImageTk
        self.canvas.delete("all")
        max_width = max(100, self.canvas.winfo_width() or self.args.canvas_width)
        max_height = max(100, self.canvas.winfo_height() or self.args.canvas_height)
        self.scale = min(max_width / self.image.width, max_height / self.image.height, 1.0)
        display = self.image.resize(
            (max(1, round(self.image.width * self.scale)), max(1, round(self.image.height * self.scale))),
            Image.Resampling.LANCZOS,
        )
        self.tk_image = ImageTk.PhotoImage(display)
        self.canvas.create_image(0, 0, anchor="nw", image=self.tk_image)
        for index, box in enumerate(self.boxes, start=1):
            self.canvas.create_rectangle(
                box.x1 * self.scale, box.y1 * self.scale,
                box.x2 * self.scale, box.y2 * self.scale,
                outline="#ff3333", width=2,
            )
            self.canvas.create_text(
                box.x1 * self.scale + 4, box.y1 * self.scale + 4,
                anchor="nw", text=str(index), fill="#ffff00",
            )

    def on_press(self, event: Any) -> None:
        if self.assessability.get() != "assessable":
            return
        self.drag_start = (event.x, event.y)
        self.drag_item = self.canvas.create_rectangle(
            event.x, event.y, event.x, event.y, outline="#00ffff", width=2
        )

    def on_drag(self, event: Any) -> None:
        if self.drag_start is not None and self.drag_item is not None:
            self.canvas.coords(self.drag_item, self.drag_start[0], self.drag_start[1], event.x, event.y)

    def on_release(self, event: Any) -> None:
        if self.drag_start is None:
            return
        x1, x2 = sorted((self.drag_start[0] / self.scale, event.x / self.scale))
        y1, y2 = sorted((self.drag_start[1] / self.scale, event.y / self.scale))
        self.drag_start = None
        self.drag_item = None
        if self.image is None:
            return
        x1, x2 = max(0, x1), min(self.image.width, x2)
        y1, y2 = max(0, y1), min(self.image.height, y2)
        if x2 - x1 < 3 or y2 - y1 < 3:
            self.render()
            return
        self.boxes.append(
            BoxRecord(
                x1, y1, x2, y2, self.life_stage.get(), self.occlusion.get(),
                bool(self.edge_truncated.get()),
            )
        )
        self.edge_truncated.set(False)
        self.render()

    def undo_box(self) -> None:
        if self.boxes:
            self.boxes.pop()
            self.render()

    def clear_boxes(self) -> None:
        self.boxes = []
        self.render()

    def rows_for_current(self) -> list[dict[str, Any]]:
        if self.image is None:
            raise RuntimeError("No image loaded")
        audit_id = text(self.tasks.iloc[self.index]["audit_id"])
        state = self.assessability.get()
        if state == "assessable" and not self.boxes:
            raise ValueError("Assessable images require at least one visible-capitulum box")
        if state != "assessable" and self.boxes:
            raise ValueError(f"{state} images cannot retain boxes")
        common = {
            "audit_id": audit_id,
            "annotator_id": self.annotator_id,
            "annotation_status": "complete",
            "assessability": state,
            "image_quality": self.image_quality.get(),
            "image_width": self.image.width,
            "image_height": self.image.height,
            "notes": self.notes.get(),
        }
        if not self.boxes:
            return [{
                **common, "gt_index": "", "gt_label": "", "life_stage": "",
                "occlusion": "", "edge_truncated": "", "gt_x1": "", "gt_y1": "",
                "gt_x2": "", "gt_y2": "",
            }]
        return [
            {
                **common,
                "gt_index": index,
                "gt_label": "visible_capitulum",
                "life_stage": box.life_stage,
                "occlusion": box.occlusion,
                "edge_truncated": box.edge_truncated,
                "gt_x1": round(box.x1, 3), "gt_y1": round(box.y1, 3),
                "gt_x2": round(box.x2, 3), "gt_y2": round(box.y2, 3),
            }
            for index, box in enumerate(self.boxes, start=1)
        ]

    def atomic_write(self) -> None:
        self.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.output.with_suffix(self.output.suffix + ".tmp")
        with temporary.open("w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(self.saved_rows)
        os.replace(temporary, self.output)

    def save_current(self) -> bool:
        audit_id = text(self.tasks.iloc[self.index]["audit_id"])
        try:
            rows = self.rows_for_current()
        except ValueError as error:
            self.messagebox.showerror("Cannot save", str(error))
            return False
        self.saved_rows = [row for row in self.saved_rows if text(row.get("audit_id")) != audit_id]
        self.saved_rows.extend(rows)
        self.saved_by_id = self._group_saved(self.saved_rows)
        self.atomic_write()
        self.progress_label.configure(
            text=f"{self.index + 1}/{len(self.tasks)} — {audit_id} — saved {len(self.saved_by_id)}/{len(self.tasks)}"
        )
        return True

    def save_next(self) -> None:
        if not self.save_current():
            return
        if self.index < len(self.tasks) - 1:
            self.index += 1
            self.load_task()
        else:
            self.messagebox.showinfo("Complete", f"All assigned images are saved to {self.output}")

    def previous(self) -> None:
        if self.index > 0:
            self.index -= 1
            self.load_task()

    def on_close(self) -> None:
        if self.messagebox.askyesno("Exit", "Exit the annotation app? Saved tasks are preserved."):
            self.root.destroy()


def main() -> None:
    args = parse_args()
    import tkinter as tk
    root = tk.Tk()
    AnnotationApp(root, args)
    root.mainloop()


if __name__ == "__main__":
    main()
