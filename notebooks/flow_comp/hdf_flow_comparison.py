import hdf_flow_comparison_functions as hf

## User inputs and working directories
print("\n\n")
inputs_filepath = fr"{input('Enter folder path with inputs:')}"

def main():
    # Get HUC10 from folder name
    huc10 = hf.get_huc_name(inputs_filepath)
    print(f"Reading HUC10 {huc10}...")

    # Create output folder if it doesn't exist
    outputs_filepath = hf.create_output_dir(inputs_filepath)

    # Get list of hdf files
    hdf_files_list = hf.get_hdf_files_with_ext(inputs_filepath)

    # List the hdf files with existing reference lines
    plans_w_ref_list = hf.check_for_ref_lines_in_hdfs(hdf_files_list, outputs_filepath)

    # Generate reference vs junction line plots for each plan and create summary table for entire model
    hf.create_plots_and_table(inputs_filepath, outputs_filepath, huc10, plans_w_ref_list)
    print ("\n\nDONE\n\n")
    return

if __name__ == "__main__":
    main()