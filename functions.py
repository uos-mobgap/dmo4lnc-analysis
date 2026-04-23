import pandas as pd
import os
from datetime import datetime
from scipy.io import savemat

from mobgap.data import GenericMobilisedDataset


def get_day_boundaries(file_path):
    """Return (start, end) row indices for each day in the CSV."""
    time_col = pd.read_csv(file_path, usecols=["Time"])
    time_col["Time"] = pd.to_datetime(time_col["Time"], format="%Y-%m-%d %H:%M:%S.%f")

    days = time_col["Time"].dt.date
    change_idx = days.ne(days.shift()).to_numpy().nonzero()[0]
    change_idx = list(change_idx) + [len(days)]

    return [(change_idx[i], change_idx[i + 1]) for i in range(len(change_idx) - 1)]


def extract_day(file_path, start, end):
    """Read only rows between [start, end) from CSV."""
    nrows = end - start
    skip = (j for j in range(1, start + 1))  # generator → memory efficient
    df = pd.read_csv(
        file_path,
        skiprows=skip,
        nrows=nrows,
        parse_dates=['Time'],  # Replace with your actual column name/index
        date_format="%Y-%m-%d %H:%M:%S.%f"
    )
    return df


def build_sensor_struct(df, fs, name):
    """Convert one day's dataframe into MATLAB struct for a sensor."""
    reformatted_time = (df["Time"].astype("int64") / 1e9).values.reshape(-1, 1)
    return {
        name: {
            "Fs": {"Acc": fs, "Gyr": fs},
            "Acc": df[["Accel-X (g)", " Accel-Y (g)", " Accel-Z (g)"]].values,
            "Gyr": df[[" Gyro-X (d/s)", " Gyro-Y (d/s)", " Gyro-Z (d/s)"]].values,
            "Timestamp": reformatted_time,
        }
    }


def load_metadata(raw_data_folder):
    """Read subject metadata once, return as infoForAlgo dict."""
    metadata = pd.read_excel(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(raw_data_folder))))
        + r"\subject_metadata.xlsx"
    )
    subject_id = int(os.path.basename(os.path.dirname(os.path.dirname(raw_data_folder))))
    subject = metadata.loc[metadata["Subject_ID"] == subject_id]

    return {
        "TimeMeasure1": {
            "Subject_ID": subject["Subject_ID"].values[0],
            "Cohort": subject["Cohort"].values[0],
            "Gender": subject["Gender"].values[0],
            "Handedness": subject["Handedness"].values[0],
            "Age": subject["Age"].values[0],
            "Weight": subject["Weight"].values[0],
            "Height": subject["Height"].values[0],
            "SensorHeight": subject["SensorHeight"].values[0],
            "WalkingAid_01": subject["WalkingAid_01"].values[0],
            "WalkingAid_Side": subject["WalkingAid_Side"].values[0],
            "WalkingAid_Description": subject["WalkingAid_Description"].values[0],
            "SensorType_SU": subject["SensorType_SU"].values[0],
            "SensorAttachment_SU": subject["SensorAttachment_SU"].values[0],
        }
    }


def save_to_day_mats(raw_data_folder, regenerate, prints):
    """
    Split wrist & waist CSVs into per-day data.mat files.
    Each file contains both sensors (if present), plus metadata.
    Skips days where data.mat already exists.
    """

    flag_file_path = os.path.join(raw_data_folder, "ALL_DAYS_DONE_FLAG")

    if regenerate:
        if os.path.exists(flag_file_path):
            os.remove(flag_file_path)

    if not os.path.isfile(flag_file_path):
        wrist_file = os.path.join(raw_data_folder, "wrist.resampled.csv")
        waist_file = os.path.join(raw_data_folder, "lower_back.resampled.csv")

        wrist_bounds = get_day_boundaries(wrist_file) if os.path.exists(wrist_file) else []
        waist_bounds = get_day_boundaries(waist_file) if os.path.exists(waist_file) else []

        n_days = max(len(wrist_bounds), len(waist_bounds))
        if prints:
            print(f"Detected {len(wrist_bounds)} wrist days, {len(waist_bounds)} waist days")
            print(f"Total output days: {n_days}")

        infoForAlgo = load_metadata(raw_data_folder)

        for day_idx in range(n_days):
            day_folder = os.path.join(raw_data_folder, f"Day {day_idx + 1}")
            data_path = os.path.join(day_folder, "data.mat")

            # Skip if already exists
            if os.path.exists(data_path) and not (regenerate):
                print(f"⏩ Skipping Day {day_idx + 1}, already exists.")
                continue

            sensors = {}

            # Wrist for this day (if available)
            if day_idx < len(wrist_bounds):
                start, end = wrist_bounds[day_idx]
                df = extract_day(wrist_file, start, end)
                fs = float(round(len(df) / (df["Time"].iloc[-1] - df["Time"].iloc[0]).total_seconds()))
                sensors.update(build_sensor_struct(df, fs, "Wrist"))

            # Waist for this day (if available)
            if day_idx < len(waist_bounds):
                start, end = waist_bounds[day_idx]
                df = extract_day(waist_file, start, end)
                fs = float(round(len(df) / (df["Time"].iloc[-1] - df["Time"].iloc[0]).total_seconds()))
                sensors.update(build_sensor_struct(df, fs, "LowerBack"))

            # Skip empty days
            if not sensors:
                continue

            # Pick earliest timestamp across sensors for StartDateTime - should be midnight for both so not super important if there's a small difference
            all_times = [
                s["Timestamp"][0][0] for s in sensors.values() if len(s["Timestamp"]) > 0
            ]
            start_time = datetime.fromtimestamp(min(all_times)).strftime("%d-%b-%Y %H:%M:%S")

            data = {
                "TimeMeasure1": {
                    "Test1": {
                        "Trial1": {
                            "SU": sensors,
                            "StartDateTime": start_time,
                            "TimeZone": "Europe/UK",
                        }
                    }
                }
            }

            os.makedirs(day_folder, exist_ok=True)

            # Save combined .mat files
            savemat(data_path, {"data": data})
            savemat(os.path.join(day_folder, "infoForAlgo.mat"), {"infoForAlgo": infoForAlgo})

            print(f"✅ Saved Day {day_idx + 1} with sensors: {list(sensors.keys())}") # not AI I just like emojis

        # create blank file as flag that all days are done
        with open(flag_file_path, "w") as file:
            pass


def comp(list1, list2):
    """
    Helper function to compare two lists.
    :param list1: first list
    :param list2: second list
    :return: boolean value of if an element of one list is in the other
    """
    for val in list1:
        if val in list2:
            return True
    return False

def get_paths_with_extension(extension,
                             start_location=os.getcwd(),
                             folders_to_ignore = []):
    paths = []
    for root, dirs, files in os.walk(start_location):
        if comp(folders_to_ignore, root.split(os.sep)):
            continue

        for file in files:
            if file.lower().endswith(extension):
                path = os.path.join(root, file)
                paths.append(path)

    return paths


def get_mobilised_dataset(paths_list,
                          test_level_names=["TimeMeasure", "Test", "Trial"],
                          measurement_condition="laboratory",
                          parent_folders_as_metadata=["cohort", "subject_id", "location", "day"]):

    print("Building dataset from .mat files...")
    dataset = GenericMobilisedDataset(
        paths_list=paths_list, # full paths to data.mat files
        test_level_names=test_level_names, # go from outermost name to innermost
        measurement_condition=measurement_condition, # can be laboratory or free_living
        parent_folders_as_metadata=parent_folders_as_metadata # go from outermost folder to innermost
    )

    return dataset