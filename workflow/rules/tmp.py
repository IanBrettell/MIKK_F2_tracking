def get_final_csvs():
    # Get path of csv files
    traj_wo_gaps_file = os.path.join(
        "/hps/nobackup/birney/users/ian/MIKK_F2_tracking",
        "split/open_field/session_20211117_1326_R_q1/trajectories_wo_gaps/trajectories_wo_gaps.trajectories.csv"
        ) ;
    traj_file = os.path.join(
        "/hps/nobackup/birney/users/ian/MIKK_F2_tracking",
        "split/open_field/session_20211117_1326_R_q1/trajectories/trajectories.trajectories.csv"
        )
    # If there is no "without gaps" file, return the "trajectories" file
    if os.path.exists(traj_wo_gaps_file):
        return(traj_wo_gaps_file)
    else:
        return(traj_file)

samples_df = pd.read_csv("config/samples.csv", comment = '#')


def get_ref_loc():
    target_col = "cab_coords_" + "of_" + "q4"
    ref_loc = samples_df.loc[samples_df["sample"] == "20211201_1414_R", target_col]
#    ref_loc = ref_loc.to_string()
    return(ref_loc)