from reproducibility.render_current_figure_surface import FIGURES


def test_current_geb_figure_surface_labels_are_frozen_and_unique():
    labels = [row[0] for row in FIGURES]
    assert labels[:5] == [f"Figure {i}" for i in range(1, 6)]
    assert labels[5:] == [f"Figure S1.{i}" for i in range(1, 8)] + ["Figure S2.1"]
    assert len(labels) == len(set(labels)) == 13


def test_current_geb_figure_surface_roles_match_frozen_plan():
    by_label = {label: (stem, kind, source) for label, stem, kind, source in FIGURES}
    assert by_label["Figure 1"][1] == "layout"
    assert by_label["Figure 2"][2] == "Figure_2_v2_geographic_sampling_domain"
    assert by_label["Figure 3"][1] == "scale"
    assert by_label["Figure 4"][1] == "construct"
    assert by_label["Figure 5"][2] == "Figure_5_v2_candidate_robustness"
    assert by_label["Figure S1.5"][1] == "layout"
    assert by_label["Figure S2.1"][2] == "Figure_3_v2_taxon_mean_information_loss"
