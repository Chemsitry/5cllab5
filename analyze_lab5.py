"""Reproduce Lab 5 hydrogen analysis from the preserved main.tex notes.

The supplied Atomic Spectra handout labels itself Lab 6 internally. The
report follows the user's Lab 5 numbering and the handout's experiment
numbers. Experiment 3 is qualitative under the flood-modified assignment.

The sample SD of four helium D3 placements is used for a new line position,
and its standard error is used for their mean, as requested by the handout.
Bench position uncertainties were not recorded, so separations are held
fixed. The covariance includes the common calculated grating spacing.
No artificial repeat observations or quantitative mystery-gas data are used.
"""
from pathlib import Path
import csv
import json
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
FIGURES = ROOT / 'figures'
DATA = ROOT / 'data'


def run_analysis():
    FIGURES.mkdir(exist_ok=True)
    DATA.mkdir(exist_ok=True)
    helium_y = np.array([11.11, 11.18, 11.30, 11.30])
    mean_y = float(np.mean(helium_y))
    sd_y = float(np.std(helium_y, ddof=1))
    sem_y = sd_y / np.sqrt(len(helium_y))
    length_he = 45.71 - 9.13
    reference_he = 587.6
    spacing = reference_he * np.hypot(length_he, mean_y) / mean_y
    spacing_u = spacing * length_he**2 / (mean_y * (length_he**2 + mean_y**2)) * sem_y

    names = ['Red', 'Cyan', 'Blue/violet']
    colors = ['tab:red', 'tab:cyan', 'tab:purple']
    offsets = np.array([9.78, 6.97, 6.24])
    length_h = 42.20 - 14.35
    ni = np.array([3., 4., 5.])
    reference_rydberg = 1.0973731568508e7
    # Evaluate the handout's Balmer formula before rounding for display.
    reference_h = 1e9 / (reference_rydberg * (.25 - 1/ni**2))
    wavelength = spacing * offsets / np.hypot(length_h, offsets)
    derivative_d = wavelength / spacing
    derivative_y = spacing * length_h**2 / (length_h**2 + offsets**2)**1.5
    wave_cov = (np.outer(derivative_d, derivative_d) * spacing_u**2
                + np.diag((derivative_y * sd_y)**2))
    wave_u = np.sqrt(np.diag(wave_cov))

    # Propagate into wavenumber. Off-diagonal terms retain the shared d error.
    wavenumber = 1e9 / wavelength
    derivative_wave = -1e9 / wavelength**2
    wave_number_cov = np.outer(derivative_wave, derivative_wave) * wave_cov
    wave_number_u = np.sqrt(np.diag(wave_number_cov))
    coordinate = 1 / ni**2
    design = np.column_stack([coordinate, np.ones(3)])
    weight = np.linalg.inv(wave_number_cov)
    parameter_cov = np.linalg.inv(design.T @ weight @ design)
    slope, intercept = parameter_cov @ design.T @ weight @ wavenumber
    residual = wavenumber - (slope * coordinate + intercept)
    chi2 = float(residual @ weight @ residual)
    rydberg = -float(slope)
    rydberg_u = float(np.sqrt(parameter_cov[0, 0]))
    final_level = float(np.sqrt(-slope / intercept))
    final_gradient = np.array([-1/(2*intercept*final_level),
                               slope/(2*intercept**2*final_level)])
    final_level_u = float(np.sqrt(final_gradient @ parameter_cov @ final_gradient))
    wavelength_difference = wavelength - reference_h
    line_deviation = wavelength_difference / wave_u
    fixed_model_chi2 = float(wavelength_difference @ np.linalg.solve(wave_cov,wavelength_difference))
    # Chi-square survival functions for exactly one and three degrees of
    # freedom have these closed forms. They are conditional on treating the
    # estimated reading scales as fixed Gaussian standard uncertainties.
    fit_p = math.erfc(math.sqrt(chi2/2))
    fixed_p = (math.erfc(math.sqrt(fixed_model_chi2/2))
               + math.sqrt(2*fixed_model_chi2/math.pi)*math.exp(-fixed_model_chi2/2))

    helium_other_y = np.array([12.87, mean_y, 9.51, 8.46])
    helium_wave = spacing * helium_other_y / np.hypot(length_he, helium_other_y)
    helium_wave_u = np.hypot(helium_wave / spacing * spacing_u,
                             spacing * length_he**2 /
                             (length_he**2 + helium_other_y**2)**1.5 * sd_y)
    # The D3 wavelength defines d, so back-substitution of the same mean is
    # identically the fixed reference and does not have independent errors.
    helium_wave_u[1] = 0.
    result = {
        'helium': {'D3_offsets_cm':helium_y.tolist(), 'sample_count':4,
                   'mean_offset_cm':mean_y, 'sample_sd_cm':sd_y,
                   'mean_standard_error_cm':float(sem_y), 'distance_cm':length_he,
                   'reference_D3_nm':reference_he, 'spacing_nm':float(spacing),
                   'spacing_u_nm':float(spacing_u),
                   'other_lines_nm':helium_wave.tolist(),
                   'other_lines_u_nm':helium_wave_u.tolist()},
        'hydrogen': {'distance_cm':length_h, 'colors':names,
                     'offsets_cm':offsets.tolist(), 'offset_u_cm':sd_y,
                     'wavelength_nm':wavelength.tolist(), 'wavelength_u_nm':wave_u.tolist(),
                     'wavelength_covariance_nm2':wave_cov.tolist(),
                     'initial_levels':ni.tolist(), 'reference_wavelength_nm':reference_h.tolist(),
                     'wavelength_difference_nm':wavelength_difference.tolist(),
                     'individual_reference_deviation':line_deviation.tolist(),
                     'spacing_error_component_nm':(derivative_d*spacing_u).tolist(),
                     'position_error_component_nm':(derivative_y*sd_y).tolist(),
                     'inverse_initial_level_squared':coordinate.tolist(),
                     'wavenumber_per_m':wavenumber.tolist(),
                     'wavenumber_u_per_m':wave_number_u.tolist(),
                     'wavenumber_covariance':wave_number_cov.tolist()},
        'fit': {'model':'1/lambda = slope/ni^2 + intercept, with both parameters free',
                'slope_per_m':float(slope), 'intercept_per_m':float(intercept),
                'parameter_covariance':parameter_cov.tolist(),
                'R_per_m':rydberg, 'R_u_per_m':rydberg_u,
                'final_level':final_level, 'final_level_u':final_level_u,
                'chi_squared':chi2, 'degrees_of_freedom':1,
                'reduced_chi_squared':chi2, 'conditional_p_value':fit_p,
                'slope_u_per_m':float(np.sqrt(parameter_cov[0,0])),
                'intercept_u_per_m':float(np.sqrt(parameter_cov[1,1])),
                'raw_residual_per_m':residual.tolist(),
                'input_uncertainties_rescaled':False},
        'agreement': {'reference_R_per_m':reference_rydberg,
                      'R_deviation_in_u':(rydberg-reference_rydberg)/rydberg_u,
                      'nf_deviation_in_u':(final_level-2)/final_level_u,
                      'fixed_Balmer_chi_squared':fixed_model_chi2,
                      'fixed_Balmer_degrees_of_freedom':3,
                      'fixed_Balmer_conditional_p_value':fixed_p,
                      'p_value_assumption':'Gaussian reading errors with estimated scales held fixed'},
        'assumptions': {'bench_separations_held_fixed':True,
                        'reason':'Bench position uncertainties were not recorded',
                        'hydrogen_line_u_source':'Sample SD of the four helium D3 placements',
                        'hydrogen_mirror_trace_recorded':False},
        'experiment3': {'mode':'Qualitative modified lab after flooding',
                        'quantitative_analysis':False,
                        'B_observation':'Prominent red, orange, and yellow features',
                        'B_candidate':'Neon',
                        'C_observation':'Many lines, including blue, two green regions, and yellow to red',
                        'C_candidate':'Argon',
                        'identification_limit':'Qualitative candidates, without wavelength validation'},
        'sources': {'handout':'LabON05-AtomicSpectra-v25.1.0.pdf',
                    'hydrogen_reference':'Handout Rydberg constant and Balmer formula, evaluated before rounding',
                    'neon':'https://physics.nist.gov/PhysRefData/Handbook/Tables/neontable2.htm',
                    'argon':'https://physics.nist.gov/PhysRefData/Handbook/Tables/argontable2.htm'}}
    (ROOT/'analysis_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    with (DATA/'hydrogen.csv').open('w',newline='',encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['color','offset_cm','offset_u_cm','wavelength_nm','wavelength_u_nm',
                         'initial_n','wavenumber_per_m','wavenumber_u_per_m','raw_residual_per_m',
                         'predicted_wavelength_nm','difference_nm','individual_reference_deviation'])
        writer.writerows(zip(names,offsets,[sd_y]*3,wavelength,wave_u,ni,
                             wavenumber,wave_number_u,residual,
                             reference_h,wavelength_difference,line_deviation))
    with (DATA/'helium_D3.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream)
        writer.writerow(['observer','side','offset_cm','original_estimated_u_cm'])
        # main.tex says Isaac right 11.11, left 11.18. The partner table swaps sides.
        writer.writerows([['Isaac','Right',11.11,.15],['Isaac','Left',11.18,.15],
                          ['Forrest','Right',11.3,.1],['Forrest','Left',11.3,.1]])

    plt.rcParams.update({'font.size':10, 'axes.spines.top':False,
                         'axes.spines.right':False, 'savefig.facecolor':'white'})
    def save(fig,name):
        for extension in ['pdf','png']:
            fig.savefig(FIGURES/f'{name}.{extension}',dpi=220,bbox_inches='tight')
        plt.close(fig)

    # Separate full-width figures keep the photo, spectrum, and fit legible.
    fig,spectrum=plt.subplots(figsize=(7.2,1.55),layout='constrained')
    for i,color in enumerate(colors):
        spectrum.vlines(reference_h[i],-.10,.10,color=color,lw=2)
        spectrum.errorbar(wavelength[i],1,xerr=wave_u[i],fmt='|',ms=20,
                          capsize=4,color=color)
        spectrum.text(wavelength[i],1.23,f'{int(ni[i])} → 2',ha='center',fontsize=9)
    spectrum.set(xlim=(400,700),ylim=(-.35,1.52),yticks=[0,1],
                 yticklabels=['Balmer reference','Measured'],xlabel='Wavelength / nm')
    spectrum.tick_params(axis='y',length=0)
    spectrum.grid(axis='x',alpha=.2)
    save(fig,'hydrogen_spectrum')

    fig,(fit_ax,residual_ax)=plt.subplots(2,1,figsize=(7.2,3.0),sharex=True,
                                         layout='constrained',
                                         gridspec_kw={'height_ratios':[2,1]})
    for i,color in enumerate(colors):
        fit_ax.errorbar(coordinate[i],wavenumber[i]/1e6,yerr=wave_number_u[i]/1e6,
                        fmt='o',color=color,capsize=3,ms=4)
        residual_ax.errorbar(coordinate[i],residual[i]/1e6,yerr=wave_number_u[i]/1e6,
                             fmt='o',color=color,capsize=3,ms=4)
    x_line=np.linspace(.03,.12,150)
    fit_ax.plot(x_line,(slope*x_line+intercept)/1e6,color='0.3',lw=1.2)
    fit_ax.set(ylabel=r'$1/\lambda$ / $10^6$ m$^{-1}$')
    fit_ax.tick_params(labelbottom=False)
    residual_ax.axhline(0,color='0.4',lw=.8)
    residual_ax.set(xlabel=r'$1/n_i^2$',ylabel=r'Residual / $10^6$ m$^{-1}$',
                    xticks=[.04,.06,.08,.10,.12])
    for ax in (fit_ax,residual_ax): ax.grid(alpha=.2)
    save(fig,'hydrogen_fit')

    # Geometry of the apparent first-order line on the spectrometer scale.
    # The inclined dashed line is a sight-line construction, not a projected
    # light ray from the central source onto the ruler.
    fig,ax=plt.subplots(figsize=(7.2,1.7),layout='constrained')
    base=8.
    height=2.1
    ax.plot([0,0],[-.55,.55],lw=4,color='0.35')
    ax.plot([base,base],[-.7,2.75],lw=2,color='0.35')
    ax.plot([0,base],[0,0],lw=1,color='0.5')
    ax.plot([0,base],[0,height],'--',lw=1.4,color='tab:red')
    ax.scatter([base,base],[0,height],s=25,color=['0.3','tab:red'],zorder=4)
    ax.annotate('',xy=(base-.45,height),xytext=(base-.45,0),
                arrowprops={'arrowstyle':'<->','color':'tab:red'})
    ax.text(base-.65,height/2,'y',ha='right',va='center',color='tab:red')
    ax.annotate('',xy=(base,-.7),xytext=(0,-.7),arrowprops={'arrowstyle':'<->'})
    ax.text(base/2,-.92,r'$L_H=27.85$ cm',ha='center',va='top')
    ax.text(0,-1.15,'Grating / eye',ha='center')
    ax.text(base+.25,0,'Central opening',va='center')
    ax.text(base+.25,height,'Apparent line',va='center')
    ax.text(base,2.95,'Scale',ha='center')
    theta=np.linspace(0,math.atan2(height,base),50)
    ax.plot(1.5*np.cos(theta),1.5*np.sin(theta),color='0.4',lw=.8)
    ax.text(1.9,.38,r'$\theta$')
    ax.set(xlim=(-1.2,11),ylim=(-1.65,3.25))
    ax.axis('off')
    save(fig,'hydrogen_geometry')
    return result


if __name__ == '__main__':
    print(json.dumps(run_analysis(),indent=2))
