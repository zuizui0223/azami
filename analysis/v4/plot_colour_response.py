"""Display saved colour effects without refitting or fabricating colour swatches."""
import json
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analysis.v4.environment_response import ENV_LABELS

ENGLISH={'chelsa_bio04':'Temperature seasonality','chelsa_bio12':'Annual precipitation',
         'chelsa_bio15':'Precipitation seasonality','chelsa_rsds_mean':'Solar radiation',
         'chelsa_vpd_mean':'VPD','chelsa_sfcwind_mean':'Wind speed',
         'chelsa_gsp':'Growing-season precipitation','chelsa_npp':'NPP'}


def plot(result,output):
    model=result['models']['colour']
    variables=result['predictors']
    fig,axes=plt.subplots(1,3,figsize=(13,8),sharey=True)
    for ax in axes:
        ax.axvline(0,color='#999999',linewidth=.8)
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='x',alpha=.15)
    for i,v in enumerate(variables):
        for suffix,offset,color,label in [
            ('__between',-.14,'#245C91','Between species'),
            ('__within',.14,'#BC591D','Within species')]:
            term=v+suffix;y=i+offset
            for j,endpoint in enumerate(('corolla_lab_lightness','corolla_lab_chroma')):
                effect=model['components'][endpoint]['effects'][term]
                b=effect['raw_units_per_environment_sd'];lo,hi=effect['raw_ci95']
                axes[j].errorbar(b,y,xerr=[[b-lo],[hi-b]],fmt='o',color=color,
                                markersize=4,capsize=2,label=label if i==0 else None)
            hue=model.get('conditional_hue_contrasts',{}).get(term,{}).get('change_degrees')
            if hue is not None:
                axes[2].plot(hue,y,'o',color=color,markersize=4)
    axes[0].set_yticks(range(len(variables)),[ENGLISH.get(v,v) for v in variables])
    axes[0].invert_yaxis()
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.94),ncol=2,frameon=False)
    axes[0].set_xlabel('Darker  <-  change in L*  ->  Lighter')
    axes[1].set_xlabel('Less vivid  <-  change in C*  ->  More vivid')
    axes[2].set_xlabel('HSV hue displacement (degrees)\nDescriptive point contrast; no interval')
    axes[0].set_title('Lightness',pad=38)
    axes[1].set_title('Chroma',pad=38)
    axes[2].set_title('Hue direction',pad=38)
    fig.suptitle('Colour response per +1 native-population SD of environment',y=.98)
    fig.text(.02,.025,'All eight selected environments fitted simultaneously; no refit. '
             'Bars: pointwise 95% intervals, not whole-colour significance.\n'
             'Hue: fixed reference contrast (other environments at reference values), '
             'not a predicted colour swatch or pigment concentration.',fontsize=9)
    fig.tight_layout(rect=(0,.09,1,.87))
    output=Path(output)
    fig.savefig(output/'colour_environment_response.png',dpi=180)
    fig.savefig(output/'colour_environment_response.pdf')
    plt.close(fig)


if __name__=='__main__':
    plot(json.loads(Path(sys.argv[1]).read_text()),Path(sys.argv[2]))
