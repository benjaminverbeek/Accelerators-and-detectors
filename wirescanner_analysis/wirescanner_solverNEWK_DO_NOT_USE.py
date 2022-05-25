# 22/5 2022
# Accelerators and Detectors, Calculation of beam emmitance

# TODO: 
# - Check calculation of k2 in M_system. NOTE: NOT DONE
# - Change matrix for quadrupole (its electric, not magnetic). NOTE: DONE
# - Implement err_sigma_11_list (real errors). NOTE: DONE
# - Implement the following (NOTE: DONE)
#   After a series of measurements, it was discovered that:
#   • The voltage that is seen by the beam in quadrupole 1 and 3 is only 0.83 times the voltage reported by the control system.
#   • Quadrupole 2 is 0.84 times weaker than the value reported by the control system.

# Imports
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from scipy.signal import find_peaks

# Variables
Z = 1 # Charge of deutron
R = 0.022 # Aperture radius for quadrupole [m]
phi = 2e6 # Acceleration voltage [V]
LQ1 = 0.1537
LQ2 = 0.287
LQ3 = 0.1537
LD1 = 0.0129
LD2 = 0.0129
LD3 = 2.3936

STD_sx2 = None
STD_sy2 = None

# Loop over the 3 measurement series
for MEASUREMENT_SERIES in ["ERROR_ESTIMATION", 1, 2, 3]: 
# for MEASUREMENT_SERIES in ["ERROR_ESTIMATION", 1]: 
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    print(f"\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    print(f"Measurements series {MEASUREMENT_SERIES}")
    print(f"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")

    # Variables
    WIRESCANNER_OFFSET = -3e-2 # NOTE: MIGHT NEED TO CHANGE THIS?
    FIDUCIAL_CALIBRATION = 0.06 # From manual, fiducial calibration [m]
    voltage = 14000 # The voltage we are analyzing [V]

    # Define file names
    if MEASUREMENT_SERIES == 1:
        filenames_voltages = [ # Filename, XVC, YVC
            ("ERR02", 15000, 15000),
            ("ERR03", 15250, 15000),
            ("ERR04", 15500, 15000),
            ("ERR05", 15750, 15000),
            ("ERR06", 16000, 15000),
            ("ERR07", 16250, 15000),
            ("ERR08", 14750, 15000),
            ("ERR09", 14500, 15000),
            ("ERR10", 14250, 15000),
            ("ERR11", 14000, 15000),
            ("ERR12", 13750, 15000),
        ]
    elif MEASUREMENT_SERIES == 2:
        filenames_voltages = [ # Filename, XVC, YVC
            ("ERR02", 15000, 15000),
            ("ERR13", 15000, 15250),
            ("ERR14", 15000, 15500),
            ("ERR15", 15000, 15750),
            ("ERR16", 15000, 16000),
            ("ERR17", 15000, 16250),
            ("ERR18", 15000, 14750),
            ("ERR19", 15000, 14500),
            ("ERR20", 15000, 14250),
            ("ERR21", 15000, 14000),
            ("ERR22", 15000, 13750),
        ]
    elif MEASUREMENT_SERIES == 3:
        filenames_voltages = [ # Filename, XVC, YVC
            ("ERR23", 13750, 13750), 
            ("ERR24", 14000, 14000),
            ("ERR25", 14250, 14250),
            ("ERR26", 14500, 14500),
            ("ERR27", 14750, 14750),
            ("ERR02", 15000, 15000),
            ("ERR28", 15250, 15250),
            ("ERR29", 15500, 15500),
            ("ERR30", 15750, 15750),
            ("ERR31", 16000, 16000),
            ("ERR32", 16250, 16250),
        ]
    elif MEASUREMENT_SERIES == "ERROR_ESTIMATION":
        filenames_voltages = [ # Filename, XVC, YVC
            ("ERR33", 15000, 15000),
            ("ERR34", 15000, 15000),
            ("ERR35", 15000, 15000),
            ("ERR36", 15000, 15000),
            ("ERR37", 15000, 15000),
            ("ERR38", 15000, 15000),
            ("ERR39", 15000, 15000),
            ("ERR40", 15000, 15000),
            ("ERR41", 15000, 15000),
            ("ERR42", 15000, 15000),
        ]
    else:
        raise Exception(f"MEASUREMENT_SERIES = {MEASUREMENT_SERIES} not supported")

    results_list = [] # Initiate a list for the results (sx^2, sy^2, XVC, YVC)

    for filename_voltage in filenames_voltages:
        filename, XVC, YVC = filename_voltage

        plots_dir = f"plots/MEASUREMENT_SERIES={MEASUREMENT_SERIES}"
        plots_filename = f"wirescanner_XVC={XVC}_YVC={YVC}"

        # Define datapaths and files
        wirescanner_data_path = f"../lab-data/beam-characterization"
        wirescanner_file = f"{wirescanner_data_path}/{filename}.CSV"
        fiducials_file = f"{wirescanner_data_path}/ERR01.CSV"

        # Read file contents
        with open(wirescanner_file) as f:
            wirescanner_lines = f.readlines()

        with open(fiducials_file) as f:
            fiducials_lines = f.readlines()

        # Extract data to arrays
        t_wirescanner_list = np.array([float(line.strip().split(",")[0]) for line in wirescanner_lines[1:]])
        v_wirescanner_list = np.array([float(line.strip().split(",")[1]) for line in wirescanner_lines[1:]])

        t_fiducials_list = np.array([float(line.strip().split(",")[0]) for line in fiducials_lines[1:]])
        v_fiducials_list = np.array([float(line.strip().split(",")[1]) for line in fiducials_lines[1:]])

        fig, (ax1, ax2) = plt.subplots(2)
        ax1.plot(t_wirescanner_list, v_wirescanner_list, linewidth=1, label="data")
        ax2.plot(t_fiducials_list, v_fiducials_list, linewidth=1, label="data")
        ax1.set_title("Wire scanner data")
        ax2.set_title("Fiducial data")
        ax1.set_xlabel("t (s)")
        ax1.set_ylabel("v (V)")
        ax2.set_ylabel("v (V)")


        # Find peaks of (-) fiducal data to get the time information
        fiducials_peak_idxs, _ = find_peaks(-v_fiducials_list, height=2, distance=1000)
        if len(fiducials_peak_idxs) == 0:
            raise Exception("No fiducial peaks found. Adjust find_peaks arguments!")
        elif len(fiducials_peak_idxs) == 1:
            raise Exception("Only one fiducial peak found. Adjust find_peaks arguments!")
        elif len(fiducials_peak_idxs) != 3:
            raise Exception(f"Found {len(fiducials_peak_idxs)} fiducial peaks (expected 3). Adjust find_peaks arguments!")
        for idx in fiducials_peak_idxs:
            ax2.plot(t_fiducials_list[idx], v_fiducials_list[idx], '*', color="purple")

        # Define minimum distance between wire scanner peaks (in indecies)
        # and also calculate the spin rate
        min_peak_distance = fiducials_peak_idxs[1] - fiducials_peak_idxs[0]
        time_between_fiducials = t_fiducials_list[fiducials_peak_idxs[1]] - t_fiducials_list[fiducials_peak_idxs[0]]
        spin_rate = 1/4 * 1/(time_between_fiducials)
        print(f"Spin rate of wire scanner: {spin_rate:.2f} cps")

        # Find peaks of wire scanner data 
        wirescanner_peak_idxs, wirescanner_peak_props = find_peaks(v_wirescanner_list, distance=1.5*min_peak_distance, width=0)
        if len(wirescanner_peak_idxs) == 0:
            raise Exception("No wirescanner peaks found. Adjust find_peaks arguments!")
        elif len(wirescanner_peak_idxs) == 1:
            raise Exception("One one wirescanner peak found. Adjust find_peaks arguments!")
        elif len(wirescanner_peak_idxs) != 2:
            raise Exception(f"Found {len(fiducials_peak_idxs)} wirescanner peaks (expected 2). Adjust find_peaks arguments!")
        for idx in wirescanner_peak_idxs:
            ax1.plot(t_wirescanner_list[idx], v_wirescanner_list[idx], '*', color="purple")


        plt.tight_layout()
        plt.savefig(f"{plots_dir}/{plots_filename}_total.png")

        peak_widths = []

        for n, idx in enumerate(wirescanner_peak_idxs):
            PLOT_WIDTH = 1000
            ax1.set_xlim(t_wirescanner_list[idx-PLOT_WIDTH], t_wirescanner_list[idx+PLOT_WIDTH])
            ax2.set_xlim(t_fiducials_list[idx-PLOT_WIDTH], t_fiducials_list[idx+PLOT_WIDTH])

            l = int(wirescanner_peak_props['left_ips'][n])
            r = int(wirescanner_peak_props['right_ips'][n])
            width_idx = r - l
            left_edge = idx-width_idx
            right_edge = idx+width_idx

            ax1.plot(t_wirescanner_list[left_edge], v_wirescanner_list[left_edge], '*', color="red")
            ax1.plot(t_wirescanner_list[right_edge], v_wirescanner_list[right_edge], '*', color="red")

            def gauss(x, a, x0, sigma):
                return a*np.exp(-(x-x0)**2/(2*sigma**2)) + WIRESCANNER_OFFSET

            p0 = [1, 0, 1]
            popt, pcov = curve_fit(gauss, t_wirescanner_list[left_edge:right_edge], v_wirescanner_list[left_edge:right_edge], p0=p0, bounds=([0, -1e-1, 0], [np.inf, 1e-1, 1]))
            print(popt)


            # Extract the peak width (in time) 
            peak_width_in_time = np.abs(popt[2])

            ax1.plot(t_wirescanner_list, gauss(t_wirescanner_list, *popt), 'k-', label=f'fit, $\sigma=${peak_width_in_time:.4e}')
            ax1.get_legend().remove() if ax1.get_legend() else None
            ax1.legend()

            # Calculate peak width (in length)
            peak_width_in_length = FIDUCIAL_CALIBRATION * peak_width_in_time / time_between_fiducials
            peak_widths.append(peak_width_in_length)

            plt.tight_layout()
            plt.savefig(f"{plots_dir}/{plots_filename}_peak-{n}.png")
            ax1.lines.pop() # remove the line we just drew in order to be able to plot in future

        sx = peak_widths[0]*1e3 # [m]
        sy = peak_widths[1]*1e3 # [m]

        results_list.append((sx**2, sy**2, XVC, YVC))

        # print("Peak widths:")
        # print(f"\tX:\ts_x = {(peak_widths[0]):.4e} m \t({(peak_widths[0]*1e3):.4f} mm)")
        # print(f"\t\ts_x^2 = {(peak_widths[0]**2):.4e} m^2 \t({(peak_widths[0]**2*1e6):.4f} mm^2)")
        # print(f"\tY:\ts_y = {(peak_widths[1]):.4e} m \t({(peak_widths[1]*1e3):.4f} mm)")
        # print(f"\t\ts_y^2 = {(peak_widths[1]**2):.4e} m^2 \t({(peak_widths[1]**2*1e6):.4f} mm^2)")

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    plt.close("all")

    if MEASUREMENT_SERIES == "ERROR_ESTIMATION":
        STD_sx2 = np.std([sx2 for sx2, _, _, _ in results_list])
        STD_sy2 = np.std([sy2 for _, sy2, _, _ in results_list])

        continue # Skip rest of analysis as is not needed when doing error estimation

    for AXIS in ["X", "Y"]:
        if MEASUREMENT_SERIES == 1:
            voltage_list = np.array([voltage for _, _, voltage, _ in results_list])
            _, _, other_voltage, _ = results_list[0]
        elif MEASUREMENT_SERIES == 2:
            voltage_list = np.array([voltage for _, _, _, voltage in results_list])
            _, _, _, other_voltage = results_list[0]
        elif MEASUREMENT_SERIES == 3:
            # Now the two voltages are the same, so it doesnt matter which one we take
            voltage_list = np.array([voltage for _, _, voltage, _ in results_list])
            other_voltage = None
        else:
            raise Exception(f"MEASUREMENT_SERIES = {MEASUREMENT_SERIES} not supported")


        if AXIS == "X":
            sigma_11_list = np.array([sx2 for sx2, _, _, _ in results_list])
            err_sigma_11_list = np.array([STD_sx2 for _ in sigma_11_list])
        elif AXIS == "Y":
            sigma_11_list = np.array([sy2 for _, sy2, _, _ in results_list])
            err_sigma_11_list = np.array([STD_sy2 for _ in sigma_11_list])
        else:
            raise Exception(f"Unsupported axis: {AXIS}")

        def k_fcn(U):
            return np.sqrt(Z*U/(R**2 * phi) + 0j)

        # Drift
        def M_D(l):
            return np.array([[1, l], [0, 1]])

        # Quadrupole
        def M_Q_nonsimple(l, k):
            c = np.cos(k*l)
            s = np.sin(k*l)

            return np.array([[c, 1/k*s], [-k*s, c]])

        def M_Q(l, k):
            return M_Q_nonsimple(l, k)

        # Calculate system matricies
        def M_system(XVC, YVC):

            # # Calculate k values with the corrections that were to be taken in consideration
            # k1 = k_fcn(0.83*XVC) 
            # k2 = -k_fcn(0.84*YVC) # k_fcn(-YVC) gives complex numbers, so maybe its ok to do this?
            # k3 = k_fcn(0.83*XVC)

            # if AXIS == "Y": # If we calculate vertical, we flip the values of k
            #     k1 = -k1
            #     k2 = -k2
            #     k3 = -k3

            # MQ1 = M_Q(LQ1, k1)
            # MQ2 = M_Q(LQ2, k2)
            # MQ3 = M_Q(LQ3, k3)

            # Calculate k values with the corrections that were to be taken in consideration
            k1 = k_fcn(0.83*XVC) 
            k2 = k_fcn(-0.84*YVC) 
            k3 = k_fcn(0.83*XVC)

            if AXIS == "X":
                pass
            elif AXIS == "Y":
                k1 = -k1
                k2 = -k2
                k3 = -k3
            else:
                raise Exception(f"Unsupported axis: {AXIS}")

            MQ1 = M_Q(LQ1, k1)
            MQ2 = M_Q(LQ2, k2)
            MQ3 = M_Q(LQ3, k3)

            MD1 = M_D(LD1)
            MD2 = M_D(LD2)
            MD3 = M_D(LD3)

            M = MD3 @ MQ3 @ MD2 @ MQ2 @ MD1 @ MQ1

            # Check to see if imaginary part in M
            any_nonreal_part = False
            for row in M:
                for elem in row:
                    if np.imag(elem) != 0:
                        any_nonreal_part = True
            if any_nonreal_part:
                print(f"Found non-real part in M! Is something wrong?\n{M=}")

            return np.real(M)

        # The equations to solve is (sigmaX can also be sigmaY)
        #   sigma_11(1) = sigmaX^2(1) = M_11^2*sigma_11(0) + 2*M_11*M_12*sigma_12(0) + M_12^2*sigma_22(0)

        # Now, define
        #   sigma_22(0) = A
        #   sigma_12(0) = B
        #   sigma_11(0) = C
        # Thus, we have
        #   sigma_11(1) = sigmaX^2(1) = M_11^2*C + 2*M_11*M_12*B + M_12^2*A

        # x is voltage, depending on which series we are doing
        def fcn(req, A, B, C):
            res = np.array([])
            for x in req:
                if MEASUREMENT_SERIES == 1:
                    XVC = x
                    YVC = other_voltage
                elif MEASUREMENT_SERIES == 2:
                    XVC = other_voltage
                    YVC = x
                elif MEASUREMENT_SERIES == 3:
                    XVC = x
                    YVC = x
                else:
                    raise Exception(f"MEASUREMENT_SERIES = {MEASUREMENT_SERIES} not supported")

                M = M_system(XVC=XVC, YVC=YVC)
                M_11 = M[0,0]
                M_12 = M[0,1]

                res = np.append(res, M_11**2*C + 2*M_11*M_12*B + M_12**2*A) 

            return res

        # Fit curve and plot
        # popt, pcov = curve_fit(fcn, voltage_list, sigma_11_list, sigma=err_sigma_11_list)
        popt, pcov = curve_fit(fcn, voltage_list, sigma_11_list)
        A, B, C = popt
        x_linspace = np.linspace(min(voltage_list)-1000, max(voltage_list)+1000)
        y_fit = fcn(x_linspace, *popt)

        # Calculate epsilon, beta and alpha
        epsilon_x = np.sqrt(A*C - B**2)
        alpha_x = -B/epsilon_x
        beta_x = C/epsilon_x

        if AXIS == "X":
            axis_label = "x"
        elif AXIS == "Y":
            axis_label = "y"
        else:
            raise Exception(f"Unsupported axis: {AXIS}")

        plt.errorbar(voltage_list, sigma_11_list, err_sigma_11_list, fmt="*", label="data")
        plt.plot(x_linspace, y_fit, label=fr"fit: $\epsilon_{axis_label}=${epsilon_x:.3f}, $\alpha_{axis_label}=${alpha_x:.3f}, $\beta_{axis_label}=${beta_x:.3f}")
        plt.xlabel("U [V]")
        plt.ylabel("$\sigma^2$ [mm]")
        plt.xlim(min(x_linspace), max(x_linspace))
        plt.legend()
        plt.title(f"Measurement series {MEASUREMENT_SERIES}, AXIS: {AXIS}")
        plt.savefig(f"plots/emmitance_series={MEASUREMENT_SERIES}_AXIS={AXIS}.png", dpi=300)
        plt.close()

        print(f"Results:")
        print(f"\tepsilon_{axis_label} = {epsilon_x:.3f}, alpha_{axis_label} = {alpha_x:.3f}, beta_{axis_label} = {beta_x:.3f}")
