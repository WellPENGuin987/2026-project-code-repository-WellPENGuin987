import os

import numpy as np

INIT_DIR = "initialization_files"


def _ensure_init_dir():
    if not os.path.isdir(INIT_DIR):
        os.makedirs(INIT_DIR, exist_ok=True)


def _init_fullpath(filename):
    # If user passed a full path, keep it. Otherwise use init folder.
    if os.path.isabs(filename):
        return filename
    return os.path.join(INIT_DIR, filename)


def SaveInputToFile(Box_Params, Time_Params, Sep_Particles, filename=None):
    """Save the current initialization parameters in text format."""
    if filename is None:
        filename = "initialization_saved.txt"

    if isinstance(filename, str) and filename.strip().lower() in ["n", "no", "none"]:
        print("Skipping save as requested.")
        return

    _ensure_init_dir()
    fullpath = _init_fullpath(filename)

    with open(fullpath, "w", encoding="utf-8") as f:
        f.write("# Box parameters: Len_X,Units_X,Len_Y,Units_Y,Len_Z,Units_Z\n")
        f.write(",".join([str(Box_Params[0]), str(Box_Params[1]), str(Box_Params[2]), str(Box_Params[3]), str(Box_Params[4]), str(Box_Params[5])]) + "\n")
        f.write("# Time parameters: Time,Units_Time,dt,Units_dt,T_plt,Units_T_plt\n")
        f.write(",".join([str(Time_Params[0]), str(Time_Params[1]), str(Time_Params[2]), str(Time_Params[3]), str(Time_Params[4]), str(Time_Params[5])]) + "\n")
        f.write("# Particle lines: type,N,m,D,r,Temp,Pos_dist,Vel_dist\n")
        for p in Sep_Particles:
            f.write(",".join([str(p[0]), str(p[1]), str(p[2]), str(p[3]), str(p[5]), str(p[6]), str(p[7]), str(p[8])]) + "\n")

    print(f"Initialization values saved to {fullpath}")


def _list_saved_files():
    """List available saved initialization files."""
    saved_files = [f for f in os.listdir(INIT_DIR) if os.path.isfile(os.path.join(INIT_DIR, f))]
    if saved_files:
        print("Available saved initialization files:")
        for f in saved_files:
            print("  -", f)


def _get_file_path():
    """Get and validate file path from user input."""
    while True:
        file_path = input("Enter filename or path (default folder initialization_files), or X to cancel: \n").strip()
        if file_path.lower() in ["x", "q", "quit", "exit"]:
            raise KeyboardInterrupt("Initialization from file canceled")
        if file_path == "":
            print("Empty filename provided. Please enter a valid filename.")
            continue
        return file_path


def _resolve_file_path(file_path):
    """Resolve file path to full path, checking existence."""
    fullpath = _init_fullpath(file_path)
    if not os.path.isfile(fullpath):
        if os.path.isabs(file_path) and os.path.isfile(file_path):
            fullpath = file_path
        else:
            print(f"File does not exist at {fullpath}. Please provide a valid filename or path.")
            return None
    return fullpath


def _parse_init_file(fullpath):
    """Parse initialization file and return parameters."""
    with open(fullpath, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]

    if len(lines) < 3:
        raise ValueError("Initialization file must include box, time, and at least one particle line.")

    # Parse box parameters
    box_items = [item.strip() for item in lines[0].split(",")]
    if len(box_items) != 6:
        raise ValueError("First line must have 6 values: Len_X,Units_X,Len_Y,Units_Y,Len_Z,Units_Z")
    Box_Params = [
        np.longdouble(str(box_items[0])),
        SetUnits(box_items[1]),
        np.longdouble(str(box_items[2])),
        SetUnits(box_items[3]),
        np.longdouble(str(box_items[4])),
        SetUnits(box_items[5]),
    ]

    # Parse time parameters
    time_items = [item.strip() for item in lines[1].split(",")]
    if len(time_items) != 6:
        raise ValueError("Second line must have 6 values: Time,Units_Time,dt,Units_dt,T_plt,Units_T_plt")
    Time_Params = [
        np.longdouble(str(time_items[0])),
        SetTime(time_items[1]),
        np.longdouble(str(time_items[2])),
        SetTime(time_items[3]),
        np.longdouble(str(time_items[4])),
        SetTime(time_items[5]),
    ]

    # Parse particle parameters
    Sep_Particles = []
    for line in lines[2:]:
        part_items = [item.strip() for item in line.split(",")]
        if len(part_items) != 8:
            raise ValueError("Particle lines must have 8 values: type,N,m,D,r,Temp,Pos_dist,Vel_dist")
        p_type = part_items[0]
        p_N = int(part_items[1])
        p_m = np.longdouble(str(part_items[2]))
        p_D = int(part_items[3])
        p_r = np.longdouble(str(part_items[4]))
        p_temp = np.longdouble(str(part_items[5]))
        p_pos = setPosDistType(int(part_items[6]))
        p_vel = setVelDistType(int(part_items[7]))

        inertia_data = setMomentsOfInertia(p_D)
        Sep_Particles.append([p_type, p_N, p_m, p_D, inertia_data, p_r, p_temp, p_pos, p_vel])

    return Box_Params, Time_Params, Sep_Particles


def ChooseFile():
    """Load initialization parameters from a CSV-like text file."""
    _ensure_init_dir()

    while True:
        _list_saved_files()
        file_path = _get_file_path()
        fullpath = _resolve_file_path(file_path)
        if fullpath is None:
            continue

        try:
            return _parse_init_file(fullpath)
        except Exception as e:
            print("Failed to parse initialization file:", e)
            print("Please fix the file or choose another file.")
            continue


def SetUnits(Units):
    allowed_units_lwr = ["nm", "mim", "mm", "cm", "m"]
    if Units.lower() in allowed_units_lwr:
        Units = Units.lower()
    elif Units == "a":
        Units = "A"
    else:
        raise ValueError
    return Units


def SetTime(Time):
    allowed_units_lwr = ["ns", "mis", "ms", "s"]
    if Time.lower() in allowed_units_lwr:
        Time = Time.lower()
    else:
        raise ValueError
    return Time


def check_positive_float(parameter, values):
    try:
        values = [np.longdouble(v) for v in values]
        if any(v <= 0 for v in values):
            raise ValueError(f"{parameter} must be positive")
        return values
    except (ValueError, TypeError) as e:
        print("Invalid input: ", e)
        raise


def setMomentsOfInertia(degrees):
    degrees = int(degrees)
    if degrees == 3:
        d, MoI, vib = Moments_of_inertia_point()
    elif degrees == 4:
        d, MoI, vib = Moments_of_inertia_linear()
    elif degrees == 6:
        d, MoI, vib = Moments_of_inertia_non_linear()
    elif degrees > 6:
        d, MoI, vib = Moments_of_inertia_vibrational(degrees)
    else:
        raise ValueError("Invalid number of degrees of freedom, must be at least 3")

    return {
        "degrees": d,
        "moments_of_inertia": MoI,
        "vibrational_modes": vib,
    }


def Moments_of_inertia_point():
    return 3, [0, 0, 0], []  # for point-like particles, moments of inertia are all zero


def Moments_of_inertia_linear():
    while True:
        moments_of_inertia = input("Moment of inertia (I) in kg*m^2\t format: I_x,I_y (I_z is the minimum of I_x & I_y):\n").replace(" ", "").split(",")
        if len(moments_of_inertia) != 2:
            print("Invalid number of moments of inertia, expected 2 for symmetric top particles, got {len(moments_of_inertia)}")
        else:
            moments_of_inertia = [
                np.longdouble(str(moments_of_inertia[0])),
                np.longdouble(str(moments_of_inertia[1])),
                min(np.longdouble(str(moments_of_inertia[0])), np.longdouble(str(moments_of_inertia[1]))),
            ]  # for symmetric top particles, the 3rd moment is the min of the first two
            check_positive_float("Moments of inertia", moments_of_inertia)
            break
    return 5, moments_of_inertia, []


def Moments_of_inertia_non_linear():
    while True:
        moments_of_inertia = input("Moments of inertia (I) in kg*m^2\t format: I_x,I_y,I_z:\n").replace(" ", "").split(",")
        if len(moments_of_inertia) != 3:
            print("Invalid number of moments of inertia, expected 3 for rotational degrees of freedom, got {len(moments_of_inertia)}")
        else:
            check_positive_float("Moments of inertia", moments_of_inertia)
            break
    return 6, moments_of_inertia, []


def Moments_of_inertia_vibrational(degrees):
    while True:
        moments_of_inertia = input("Moments of inertia (I) in kg*m^2\t format: I_x,I_y,I_z:\n").replace(" ", "").split(",")
        if len(moments_of_inertia) != 3:
            print("Invalid number of moments of inertia, expected 3 for rotational degrees of freedom, got {len(moments_of_inertia)}")
        else:
            check_positive_float("Moments of inertia", moments_of_inertia)
            break
    while True:
        vibrational_modes = input().replace(" ", "").split(",")
        if len(vibrational_modes) != degrees - 6:
            print(f"Invalid number of vibrational modes, expected {degrees - 6}, got {len(vibrational_modes)}")
        else:
            check_positive_float("Vibrational modes", vibrational_modes)
            break
    return degrees, moments_of_inertia, vibrational_modes


def setPosDistType(pos_dist_type):
    allowed_pos_dist_types = [1, 2, 3]
    if pos_dist_type in allowed_pos_dist_types:
        return pos_dist_type
    else:
        raise ValueError


def setVelDistType(vel_dist_type):
    allowed_vel_dist_types = [1, 2, 3]
    if vel_dist_type in allowed_vel_dist_types:
        return vel_dist_type
    else:
        raise ValueError


def InitialiseBox():
    while True:
        Init_Box = (
            input(
                "Input parameters for the box:\n"
                "Len_X/Y/Z are the dimensions of the container\n"
                "UnitsX/Y/Z: supports A (Ångstroms), nm, mim (micrometres), mm, cm, m\n"
                "Format:\tLen_X,Units_X,Len_Y,Units_Y,Len_Z,UnitsZ\n"
            )
            .replace(" ", "")
            .split(",")
        )
        print(Init_Box)

        if Init_Box == ["F"] or Init_Box == ["f"]:
            Box_Params, _, _ = ChooseFile()
            return Box_Params
        elif len(Init_Box) != 6:
            print("Invalid number of dimensions and units")
        else:
            try:
                Len_X, Units_X, Len_Y, Units_Y, Len_Z, Units_Z = Init_Box
                Len_X = np.longdouble(str(Len_X))
                Len_Y = np.longdouble(str(Len_Y))
                Len_Z = np.longdouble(str(Len_Z))
                Units_X = SetUnits(Units_X)
                Units_Y = SetUnits(Units_Y)
                Units_Z = SetUnits(Units_Z)
                Box_Params = [Len_X, Units_X, Len_Y, Units_Y, Len_Z, Units_Z]
                return Box_Params
            except (ValueError, TypeError):
                print("Invalid inputs")


def InitialiseTime():
    while True:
        Init_Time = (
            input(
                "Input paramenters for time:\n"
                "Time is the total run time\tdt is the time step between each test\n"
                "Units_Time/ dt/ T_plt: supports ns, mis (microseconds), ms, s\n"
                "T_plt is the time between data samples for plots\n"
                "Format:\tTime,Units_Time,dt,Units_dt,T_plt,Units_T_plt\n"
            )
            .replace(" ", "")
            .split(",")
        )

        if Init_Time == ["F"] or Init_Time == ["f"]:
            _, Time_Params, _ = ChooseFile()
            return Time_Params
        elif len(Init_Time) != 6:
            print("Invalid number of dimensions and units")
        else:
            try:
                Time, Units_Time, dt, Units_dt, T_plt, Units_T_plt = Init_Time
                Time = np.longdouble(str(Time))
                dt = np.longdouble(str(dt))
                T_plt = np.longdouble(str(T_plt))
                Units_Time = SetTime(Units_Time)
                Units_dt = SetTime(Units_dt)
                Units_T_plt = SetTime(Units_T_plt)
                Time_Params = [Time, Units_Time, dt, Units_dt, T_plt, Units_T_plt]
                return Time_Params
            except (ValueError, TypeError):
                print("Invalid inputs")


def InitialiseParticles():
    Sep_Particles = []
    while True:
        Init_Particles = (
            input(
                "Input parametetrs of particles (X to escape):\n"
                "type is the type of particle, used to differentiate different types of particles in the graphs\n"
                "N is the number of particles of this type\tm is the mass of the particles (kg)\n"
                "D is how many degrees of freedom these particles have: 3 translational (minimum) + 1 or 3 rotational + 1 vibrational\n"
                "r is the radius of collision of these particles (Å)\tTemp is the initial temperature of these particles (K)\n"
                "Pos_dist is the type of distribution of the initial positions of these particles (Lattice = 1, Uniform pdf = 2, Custom gradient pdf = 3)\n"
                "Vel_dist is the type of distribution of the initial velocities of these particles (Identical RMS = 1, Uniform RMS pdf = 2, Maxwell-Boltzmann = 3)\n"
                "Format:\ttype,N,m,D,r,Temp,Pos_dist,Vel_dist\n"
            )
            .replace(" ", "")
            .split(",")
        )

        if Init_Particles == ["X"] or Init_Particles == ["x"]:
            print("Initialisation complete")
            break
        elif Init_Particles == ["F"] or Init_Particles == ["f"]:
            _, _, Sep_Particles = ChooseFile()
            return Sep_Particles
        else:
            print(f"{len(Init_Particles)} parameters inputted:\n{Init_Particles}")
        if len(Init_Particles) % 8 == 0:
            Current_Particles = []  # list to store parameters of current particle type, will be added to list of all particle parameters
            try:
                Current_Particles.append(str(Init_Particles[0]))  # particle type
                Current_Particles.append(int(Init_Particles[1]))  # number of particles of this type
                Current_Particles.append(np.longdouble(str(Init_Particles[2])))  # mass of particles of this type

                # keep degrees as numeric for simulation, and keep inertia structure as extra metadata
                inertia_data = setMomentsOfInertia(Init_Particles[3])
                Current_Particles.append(int(Init_Particles[3]))  # degrees of freedom
                Current_Particles.append(inertia_data)  # moments of inertia / vibrational modes metadata

                Current_Particles.append(np.longdouble(str(Init_Particles[4])))  # radius of collision
                Current_Particles.append(np.longdouble(str(Init_Particles[5])))  # initial temperature
                Current_Particles.append(setPosDistType(int(Init_Particles[6])))  # position distribution type
                Current_Particles.append(setVelDistType(int(Init_Particles[7])))  # velocity distribution type

                Sep_Particles.append(Current_Particles)  # add current particle parameters to list of all particle parameters

            except (ValueError, TypeError) as e:
                print("Invalid inputs: ", e)
        else:
            print("Invalid number of particle parameters")
    return Sep_Particles


def Initialise():
    # --------------initialising the parameters of the simulation--------------#
    use_file = input("Use file initialization? (Y/N): ").strip().lower()
    if use_file in ["y", "yes", "f", "file"]:
        try:
            Box_Params, Time_Params, Sep_Particles = ChooseFile()
            # When loaded from file, keep current input file as source and don't auto-save.
            return (Box_Params, Time_Params, Sep_Particles)
        except KeyboardInterrupt:
            print("File initialization canceled, switching to manual input.")

    Box_Params = InitialiseBox()
    Time_Params = InitialiseTime()
    Sep_Particles = InitialiseParticles()

    # Default save path for custom manual input, unless disabled by the user
    save_choice = input("Save this custom initialization to file? (default: Enter to save, N to skip): ").strip()
    if save_choice.lower() not in ["n", "no", "none"]:
        filename = input("Enter filename or press Enter for 'initialization_saved.txt': ").strip()
        if filename == "":
            filename = "initialization_saved.txt"
        SaveInputToFile(Box_Params, Time_Params, Sep_Particles, filename)
    else:
        print("Auto-save disabled by user.")

    return (Box_Params, Time_Params, Sep_Particles)
