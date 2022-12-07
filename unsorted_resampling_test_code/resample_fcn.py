def resample(): # implemented based on down_sample() but adapted / H
    # TODO: try to fix/find/implement spline (tricubic spline) interpolators
    # Eqtools does not seem to work.
    # See if you can fix issues with scipy
    # Try ARBTools Leiken-Marsden 3D tricubic spline interpolation
    # otherwise look for another, or worst case implement yourself
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file',
                        help='File containing the unitcell to be sampled.')
    parser.add_argument('--method', type=str, 
                        choices={'linear', 'lin', 'slinear', 'nearest', 'near', 
                                  'cubic', 'quintic', 'LMtricubic', 'eqt-cubicspline', 'gsmear'}, 
                        default='linear', 
                        help='Method used to interpolate.')
    parser.add_argument('--spacing', type=float, default=0.2,
                        help='Desired spacing.')
    parser.add_argument('--factor', type=int, default=None,
                        help='Integer factor with which to multiply the number'
                        'of points in each dimension in the input grid to get '
                        'the number of point in each dimension in the predicted grid.')

    args = parser.parse_args()

    in_path = Path(args.input_file)
    if in_path.suffix == '.xsf':
        atoms = ase.io.read(in_path, read_data=False)
        data = ase.io.read(in_path, read_data=True)
    elif in_path.suffix == '.cube':
        with open(in_path, 'r') as fp:
            cube_data = ase.io.cube.read_cube(fp)
            atoms = cube_data['atoms']
            data = cube_data['data']
    else:
        raise Exception('Unknown input file type.')

    # Does data contain nans or infs?
    print('Contains NaNs: ', np.isnan(data).any())
    print('Contains Infs: ', np.isinf(data).any())
    print('Finite entries: ', np.isfinite(data).all())


    #Handle methods:
    scipy_methods = {'linear', 'lin', 'slinear', 'nearest', 'near', 
                                  'cubic', 'quintic'} 

    if args.method == 'lin' or args.method == 'linear':
        method = 'linear'
    elif args.method == 'near' or args.method == 'nearest':
        method = 'nearest'
    elif args.method == 'slinear':
        method = 'slinear'
    elif args.method == 'cubic':
        method = 'cubic'
    elif args.method == 'eqt-cubicspline':
        method = 'eqt-cubicspline'
        #boundary = 'clamped' # Boundary options are 'clamped' and 'normal' in eqtools trispline
        # Natural - sort of "linear" extrapolation, while Clamped is just repeating the outer layer once more
        boundary = 'natural'
        raise Exception("Method not available as it does not work properly.")
    elif args.method == 'LMtricubic':
        method = args.method
    elif args.method == 'quintic':
        method = 'quintic'
    elif args.method == 'gsmear':
        method = 'gsmear'
    else:
        raise Exception('Unsupported interpolation method: '+ str(args.method))


    if args.factor:
        nx, ny, nz = data.shape
        n_frac_a,  n_frac_b,  n_frac_c  = args.factor*nx, args.factor*ny, args.factor*nz
        out_path = in_path.stem + '_resampled_factor_' + str(args.factor) + '_' + str(method) + '.cube'

    else:
        n_frac_a = int(np.linalg.norm(atoms.get_cell()[:][0]) / args.spacing)
        n_frac_b = int(np.linalg.norm(atoms.get_cell()[:][1]) / args.spacing)
        n_frac_c = int(np.linalg.norm(atoms.get_cell()[:][2]) / args.spacing)
        out_path = in_path.stem + '_resampled_spacing_' + str(args.spacing) +'_' + str(method) + '.cube'

    #Maybe TODO: Lägg till periodicitet i interpolationen? "Padda" gridden

    nx, ny, nz = data.shape

    xgrid = np.linspace(0, 1, nx)
    ygrid = np.linspace(0, 1, ny)
    zgrid = np.linspace(0, 1, nz)
 
    E_thresh = 1.0e+4 # Set the upper energy threshold to use if filtering or masking
    filter_data = False
    mask_data = False
    fill_with_lower_E = True
    print('Energy threshold for filtering/masking: ', E_thresh)
    print('Filtering: ', filter_data)
    print('Masking: ', mask_data)
    print('Replace extremes with lower E value: ', fill_with_lower_E)
    if (filter_data and mask_data) or (filter_data and fill_with_lower_E) or (mask_data and fill_with_lower_E):
        raise Exception("Filtering, masking and filling with lower E values are mutually exclusive!")

    if filter_data:
        # Try filtering the data again: (see if it works)
        reduced_data = np.copy(data)
        max_data = np.max(data)
        print('Max E: ', max_data)
        reduced_data[data > E_thresh] = E_thresh
        data = reduced_data
        print('Max E after reduction: ', np.max(data))
        #
    elif mask_data:
        # Try masking data instead:
        mask = data < E_thresh
        filtered_data = np.copy(data[mask]).reshape(-1,1)
        filtered_coords = np.array(list(it.product(xgrid, ygrid, zgrid))).reshape(data.shape + (3,))[mask]
        filtered_coords = filtered_coords.reshape(-1,3)
        filtered_field = np.hstack((filtered_coords, filtered_data))
        #
    elif fill_with_lower_E:
        max_value = np.nan_as_num(np.inf)
        max_data = data == max_value
        not_max_data = data != max_value
        norm_max_val = np.max(data[not_max_data])
        data[max_data] = 10*norm_max_value

    # Construct desired interpolant
    if method != 'gsmear' and (args.method in scipy_methods):
        interpolating_function = RegularGridInterpolator((xgrid, ygrid, zgrid), data, method=method)
    elif method == 'gsmear':
        nx, ny, nz = data.shape
        #spacing = np.linalg.norm(atoms.get_cell([]))
        #half_width = 
        spaced_out_data = np.zeros((2*nx, 2*ny, 2*nz))
        spaced_out_data[0::2, 0::2, 0::2] = np.copy(data)
        data_new = scipy.ndimage.gaussian_filter(spaced_out_data, sigma=0.75, 
                order=0, output=None, mode='wrap', cval=0.0, truncate=1.5)
    elif method == 'eqt-cubicspline':
        interpolating_function = Spline(xgrid, ygrid, zgrid, data, boundary=boundary)
    elif method == 'LMtricubic':
        coord_array =  np.array(list(it.product(xgrid, ygrid, zgrid)))
        input_field = np.hstack((coord_array, data.reshape(-1,1)))
        interpolating_function = tricubic(filtered_field) 
    
    # Generate new gridpoints for prediction
    xgrid_new = np.linspace(0, 1, n_frac_a)
    ygrid_new = np.linspace(0, 1, n_frac_b)
    zgrid_new = np.linspace(0, 1, n_frac_c)
    
    X, Y, Z = np.meshgrid(xgrid_new, ygrid_new, zgrid_new)
    coord_list =  np.hstack((X.reshape(-1,1), Y.reshape(-1,1), Z.reshape(-1,1)))
    alt_coord_list = np.array(list(it.product(xgrid_new, ygrid_new, zgrid_new)))
    
    

    if method in scipy_methods:
        data_new = interpolating_function((X, Y, Z))
        #data_new = interpolating_function(coord_list)
        print("Data new shape:", data_new.shape)

    elif method == 'eqt-cubicspline':
        data_new = interpolating_function(xgrid_new, ygrid_new, zgrid_new)
        print("Data new shape:", data_new.shape)

    elif method == 'LMtricubic':
        data_new, grad = interpolating_function.Query(alt_coord_list)
        print("Data length, shape:", len(data_new), data_new.shape)
#        np.nan_to_num(data_new, copy=False)
    print(type(data_new))


    with open(out_path, 'w') as fp:
        ase.io.cube.write_cube(fp, atoms, data_new.reshape(n_frac_a, n_frac_b, n_frac_c))

