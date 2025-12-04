##

def reshape_var(da):
    if 'ncol' not in da.dims:
        return da
    # reshape preserving other dims
    dims = tuple(d for d in da.dims if d != 'ncol') + ('ny', 'nx')
    return da.values.reshape(*(da.sizes[d] for d in da.dims if d != 'ncol'), ny, nx)


def reshape_dataset(ds):
    
    ny, nx = len(ds.lat_R), len(ds.lon_R)
    # apply to all DataArrays
    for v in ds.data_vars:
        if 'ncol' in ds[v].dims:
            da = ds[v]
            new_dims = tuple(d for d in da.dims if d != 'ncol') + ('ny', 'nx')
            ds[v] = (new_dims, da.values.reshape(*(da.sizes[d] for d in da.dims if d != 'ncol'), ny, nx))
    
    # add coordinates
    ds = ds.assign_coords(lat=('ny', ds.lat_R), lon=('nx', ds.lon_R))

    return ds