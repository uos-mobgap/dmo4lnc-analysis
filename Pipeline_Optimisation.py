import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from functions import *
from mobgap.pipeline import MobilisedPipelineImpaired
from mobgap.aggregation import get_mobilised_dmo_thresholds
from mobgap.gait_sequences import GsdAdaptiveIonescu, GsdIluz, GsdIonescu
from mobgap.initial_contacts import IcdHKLeeImproved, IcdIonescu, IcdShinImproved

from skopt import gp_minimize
from skopt.space import Real, Categorical
from skopt.plots import plot_convergence, plot_evaluations, plot_objective

warnings.filterwarnings("ignore", category=UserWarning)


def calculate_loso_cost(metrics, dataset, cohorts, subjects_to_ignore):
    """Calculates mean cost function across all LOSO folds"""
    if metrics is None or metrics.empty:
        return 1e4

    # get subject list for the given cohorts
    active_subjects = []
    for cohort in cohorts:
        cohort_subs = {row[1] for row in dataset.group_labels if row[0] == cohort}
        active_subjects.extend(list(cohort_subs - set(subjects_to_ignore)))

    # if we don't have any subjects (likely incorrect cohort) return big number
    if not active_subjects:
        return 1e4

    fold_costs = []

    # evaluate cost for each fold
    for sub in active_subjects:
        sub_metrics = metrics[metrics["subject"] == sub]

        if len(sub_metrics) > 0 and sub_metrics["wb_TPs"].sum() > 0:
            wb_FPs = sub_metrics["wb_FPs"].sum()
            wb_FNs = sub_metrics["wb_FNs"].sum()

            start_time_diff = (sub_metrics["wb_start_mocap"] - sub_metrics["wb_start_mobgap"]).abs().sum()
            end_time_diff = (sub_metrics["wb_end_mocap"] - sub_metrics["wb_end_mobgap"]).abs().sum()
            duration_diff = (sub_metrics["wb_duration_mocap"] - sub_metrics["wb_duration_mobgap"]).abs().sum()

            missing_bout_weight = 10
            fold_cost = (missing_bout_weight * wb_FPs +
                         missing_bout_weight * wb_FNs +
                         start_time_diff +
                         end_time_diff +
                         duration_diff)
        else:
            fold_cost = 1e4

        fold_costs.append(fold_cost)

    return float(np.mean(fold_costs))


def assess_default_gsd():
    pipeline = MobilisedPipelineImpaired(
        gait_sequence_detection=GsdIluz(),
        dmo_thresholds=multi_cohort_thresholds
    )

    metrics = assess_pipeline(
        dataset,
        pipeline,
        cohorts,
        subjects_to_ignore,
        all_tests,
        prints=False,
    )

    return calculate_loso_cost(metrics, dataset, cohorts, subjects_to_ignore)


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
        gait_sequence_detection=GsdIluz(
            window_length_s=window_length_s,
            window_overlap=window_overlap,
            std_activity_threshold=std_activity_threshold,
            mean_activity_threshold=mean_activity_threshold,
            acc_v_standing_threshold=acc_v_standing_threshold,
            sin_template_freq_hz=sin_template_freq_hz,
            allowed_acc_v_change_per_window=allowed_acc_v_change_per_window,
            min_gsd_duration_s=min_gsd_duration_s,
            use_original_peak_detection=use_original_peak_detection
        ),
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

    mean_loso_cost = calculate_loso_cost(metrics, dataset, cohorts, subjects_to_ignore)
    print(f"Mean LOSO CV cost function = {mean_loso_cost}\n")
    return mean_loso_cost


# Get mobgap dataset
paths_list = get_paths_with_extension(
    "data.mat",
    start_location=os.path.join(os.getcwd(), "Dataset"),
    folders_to_ignore=["Home"]
)

dataset = get_mobilised_dataset(
    paths_list,
    parent_folders_as_metadata=["cohort", "subject_id", "location", "sensor"]
)
print(dataset)

cohorts = ["CP"]
subjects_to_ignore = ['969', '921']  # Excluded subjects
all_tests = ["Test" + str(i) for i in range(1, 10)]

# hack CP cohort thresholds as healthy
ha_thresholds = get_mobilised_dmo_thresholds().xs("HA", level=1, drop_level=False)
new_index = pd.MultiIndex.from_tuples(
    [(dmo, cohort) for dmo, _ in ha_thresholds.index for cohort in cohorts],
    names=ha_thresholds.index.names
)
duplicated_values = pd.concat([ha_thresholds] * len(cohorts), axis=0).reset_index(drop=True)
multi_cohort_thresholds = pd.DataFrame(duplicated_values.values, index=new_index, columns=ha_thresholds.columns)

# default cost for comparison
default_cost_function = assess_default_gsd()
print(f"Default mean LOSO cost function = {default_cost_function}\n")

space = [
    Real(0.5, 5.0, name='window_length_s'),
    Real(0.1, 0.9, name='window_overlap'),
    Real(-5.0, 5.0, name='std_activity_threshold'),
    Real(-1.0, 1.0, name='mean_activity_threshold'),
    Real(-10, 10, name='acc_v_standing_threshold'),
    Real(1, 10, name='sin_template_freq_hz'),
    Real(0.0000001, 2.0, name='allowed_acc_v_change_per_window'),
    Real(1.0, 10.0, name='min_gsd_duration_s'),
    Categorical([True, False], name='use_original_peak_detection')
]

# set to 120-200 iterations for real run
number_of_iterations = 50
number_of_random_starts = max(1, number_of_iterations // 5)

result = gp_minimize(
    gsd_eval,
    space,
    n_calls=number_of_iterations,
    n_random_starts=number_of_random_starts,
    random_state=42,
    verbose=True
)

best_params = result.x
best_score = result.fun

print(f"Best parameters: {best_params}")
print(f"Default mean LOSO cost function = {default_cost_function}")
print(f"Best minimum LOSO cost function: {best_score}")

# plot results
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
fig.subplots_adjust(hspace=0.6, wspace=0.6)
plt.show()

plot_objective(result)
fig = plt.gcf()
fig.set_size_inches(22, 22)
fig.subplots_adjust(hspace=0.6, wspace=0.6)
plt.show()

# save results
param_names = [
    dim.name if dim.name else f"param_{i}"
    for i, dim in enumerate(result.space.dimensions)
]

df_results = pd.DataFrame(result.x_iters, columns=param_names)
df_results["objective_val"] = result.func_vals

output_path = r"C:\Users\ac4jmi\Desktop\DMO4LNC\dmo4lnc-analysis\In-Lab Parameters\Pipeline Experiments\gp_minimize_results.csv"
os.makedirs(os.path.dirname(output_path), exist_ok=True)
df_results.to_csv(output_path, index=False)