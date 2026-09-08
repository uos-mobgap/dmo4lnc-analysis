from functions import *

import matplotlib.pyplot as plt

from mobgap.pipeline import MobilisedPipelineImpaired

from mobgap.aggregation import get_mobilised_dmo_thresholds
from mobgap.gait_sequences import GsdAdaptiveIonescu, GsdIluz, GsdIonescu
from mobgap.initial_contacts import IcdHKLeeImproved, IcdIonescu, IcdShinImproved

from skopt import gp_minimize
from skopt.space import Real, Categorical
from skopt.plots import plot_convergence, plot_evaluations, plot_objective

# turn off userwarnings
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

# get mobgap dataset
paths_list = get_paths_with_extension("data.mat",
                                      start_location=os.path.join(os.getcwd(), "Dataset"),
                                      folders_to_ignore=["Home"])

dataset = get_mobilised_dataset(paths_list, parent_folders_as_metadata=["cohort", "subject_id", "location", "sensor"])
print(dataset)

cohorts = ["CP"] # can add "HA" or "PSP". Optimisation should be cohort specific

subjects_to_ignore = ['969', '921'] # 969 was heavily gait impaired, 921 was not actually cp

all_tests = ["Test" + str(i) for i in range(1, 10)]

# get the HA healthy thresholds and rename them to "CP" so that we can use them to do some basic thresholding on the data
ha_thresholds = get_mobilised_dmo_thresholds().xs("HA", level=1, drop_level=False)
new_index = pd.MultiIndex.from_tuples([(dmo, cohort) for dmo, _ in ha_thresholds.index for cohort in cohorts], names=ha_thresholds.index.names)
duplicated_values = pd.concat([ha_thresholds] * len(cohorts), axis=0).reset_index(drop=True)
multi_cohort_thresholds = pd.DataFrame(duplicated_values.values, index=new_index, columns=ha_thresholds.columns)


def gsd_eval(params):
    (window_length_s, window_overlap, std_activity_threshold,
     mean_activity_threshold, acc_v_standing_threshold, sin_template_freq_hz,
     allowed_acc_v_change_per_window, min_gsd_duration_s, use_original_peak_detection) = params

    print(f"window_length_s = {window_length_s}")
    print(f"window_overlap = {window_overlap}")
    print(f"std_activity_threshold = {std_activity_threshold}")
    print(f"mean_activity_threshold = {mean_activity_threshold}")
    print(f"acc_v_standing_threshold = {acc_v_standing_threshold}")
    print(f"sin_template_freq_hz = {sin_template_freq_hz}")
    print(f"allowed_acc_v_change_per_window = {allowed_acc_v_change_per_window}")
    print(f"min_gsd_duration_s = {min_gsd_duration_s}")
    print(f"use_original_peak_detection = {use_original_peak_detection}")

    pipeline = MobilisedPipelineImpaired(
        gait_sequence_detection=GsdIluz(window_length_s = window_length_s, #
                                        window_overlap = window_overlap,
                                        std_activity_threshold = std_activity_threshold,
                                        mean_activity_threshold = mean_activity_threshold,
                                        acc_v_standing_threshold = acc_v_standing_threshold,
                                        sin_template_freq_hz = sin_template_freq_hz,
                                        allowed_acc_v_change_per_window = allowed_acc_v_change_per_window,
                                        min_gsd_duration_s = min_gsd_duration_s,
                                        use_original_peak_detection = use_original_peak_detection),
        dmo_thresholds=multi_cohort_thresholds,
    )

    metrics = assess_pipeline(
        dataset,
        pipeline,
        cohorts,
        subjects_to_ignore,
        all_tests,
        prints=False,
    )

    if metrics["wb_TPs"].sum() > 0: # if any walking bouts were detected
        wb_FPs = metrics["wb_FPs"].sum()
        wb_FNs = metrics["wb_FNs"].sum()

        start_time_diff = (metrics["wb_start_mocap"] - metrics["wb_start_mobgap"]).abs().sum()
        end_time_diff = (metrics["wb_end_mocap"] - metrics["wb_end_mobgap"]).abs().sum()
        duration_diff = (metrics["wb_duration_mocap"] - metrics["wb_end_mobgap"]).abs().sum()

        # gp_minimize minimises this value
        missing_bout_weight = 10
        cost_function = (missing_bout_weight*wb_FPs + missing_bout_weight*wb_FNs +
                         start_time_diff + end_time_diff + duration_diff)
    else:
        cost_function = 1e4
    print(f"cost function = {cost_function}")
    print()
    return cost_function

space = [
    Real(0.5, 5.0, name='window_length_s'),
    Real(0.1, 0.9, name='window_overlap'),
    Real(-5.0, 5.0, name='std_activity_threshold'),
    Real(-1.0, 1.0, name='mean_activity_threshold'),
    Real(-10, 10, name='acc_v_standing_threshold'),
    Real(1, 10, name = 'sin_template_freq_hz'),
    Real(0.0000001, 2.0, name = 'allowed_acc_v_change_per_window'),
    Real(1.0, 10.0, name = 'min_gsd_duration_s'),
    Categorical([True, False], name = 'use_original_peak_detection')
]

number_of_iterations = 120
number_of_random_starts = number_of_iterations // 5 # at 10%, this didn't explore one of the parameters properly

result = gp_minimize(
    gsd_eval,
    space,
    n_calls=number_of_iterations,
    n_random_starts=number_of_random_starts,
    random_state=42
)

best_params = result.x
best_score = result.fun

print(f"Best parameters: {best_params}")
print(f"Best minimum value: {best_score}")

# 1. Set global micro-fonts prior to plotting
plt.rcParams.update({
    'font.size': 5,
    'axes.labelsize': 5,
    'axes.titlesize': 6,
    'xtick.labelsize': 4,
    'ytick.labelsize': 4
})

plot_evaluations(result)
fig = plt.gcf()
fig.set_size_inches(22, 22)
fig.subplots_adjust(hspace=0.6, wspace=0.6)  # Increases spacing between grid items
plt.show()

plot_objective(result)
fig = plt.gcf()
fig.set_size_inches(22, 22)
fig.subplots_adjust(hspace=0.6, wspace=0.6)
plt.show()

# Save to CSV
param_names = [
    dim.name if dim.name else f"param_{i}"
    for i, dim in enumerate(result.space.dimensions)
]

df_results = pd.DataFrame(result.x_iters, columns=param_names)
df_results["objective_val"] = result.func_vals

df_results.to_csv(r"C:\Users\ac4jmi\Desktop\DMO4LNC\dmo4lnc-analysis\In-Lab Parameters\Pipeline Experiments\gp_minimize_results.csv", index=False)
