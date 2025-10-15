import numpy as np


def axes_def(n,nxplo=1,nyplo=1,xpan=0.9,ypan=0.9,xcor=1.,ycor=1.,top_panel=0,left_panel=0,top_whitespace=0.):
    dxplo=xpan/nxplo 
    dyplo=ypan/nyplo
    dxcor=xcor/nxplo
    dycor=ycor/nyplo

    jplo = (n-1)//nxplo #+ 1
    iplo = n-1 - jplo*nxplo

    top_edge = 1.0 - top_whitespace
    Axes = [ iplo*dxcor , top_edge-(jplo+1)*dycor , dxplo , dyplo ] 
    
    if ((top_panel>0)and(left_panel==0)):
        dyplo1=dyplo-top_panel
        Axes0 = [ iplo*dxcor , top_edge-(jplo+1)*dycor , dxplo , dyplo1 ] 
        Axes1 = [ iplo*dxcor , top_edge-(jplo+1)*dycor + dyplo1, dxplo , top_panel ] 
        Axes  = (Axes0,Axes1)
    if ((top_panel>0)and(left_panel>0)):
        dyplo1=dyplo-top_panel
        dxplo1=dxplo-left_panel
        Axes0 = [ iplo*dxcor        , top_edge-(jplo+1)*dycor          , dxplo1     , dyplo1 ] 
        Axes1 = [ iplo*dxcor        , top_edge-(jplo+1)*dycor + dyplo1 , dxplo1     , top_panel ] 
        Axes2 = [ iplo*dxcor+dxplo1 , top_edge-(jplo+1)*dycor          , left_panel , dyplo1 ] 
        Axes  = (Axes0,Axes1,Axes2)
    
    return Axes

def param_text_box( ax, parameters=None, xpo=0.65, ypo=0.98, select=None ):

    if (select is not None) and (parameters is not None):
        # Parameters to display in the inset box
        selected_keys = select #['alpha_gw_movmtn', 'effgw_rdg_resid','use_gw_rdg_resid']
        selected_params = {key: parameters[key] for key in selected_keys if key in parameters}
        # Format selected parameters into a string
        param_text = "\n".join([f"{key}: {value}" for key, value in selected_params.items()])
        # Add inset box
        ax.text(
            xpo, ypo,  # X, Y position in axes coordinates
            param_text,  # Text content
            transform=ax.transAxes,  # Use axes coordinates
            fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='lightgrey')
            )
    

def right_cbar_axes(Axes,xshift=0,yshift=0, xsize=0.01,ysize=0.3 ):
    CrAxes = [ Axes[0]+Axes[2]+xshift ,Axes[1]+yshift , xsize , ysize  ] 
    return CrAxes

def under_cbar_axes(Axes,xshift=0,yshift=0, xsize=0.3,ysize=0.01 ):
    CrAxes = [ Axes[0]+xshift ,Axes[1]+yshift , xsize , ysize  ] 
    return CrAxes

def step_cmap(cmapN,levels):
    import matplotlib.colors as colors
    import matplotlib.pyplot as plt


    # Normalized contour levels
    norm_levels = (levels - levels.min()) / (levels.max() - levels.min())
    # Creating a colormap from the "viridis" colormap
    colomap = plt.get_cmap(cmapN)
    colors_array = colomap(np.linspace(0, 1, len(norm_levels)))
    # Creating a new colormap
    cmap = colors.LinearSegmentedColormap.from_list("nonuniform", list(zip(norm_levels, colors_array)))

    return cmap

def save_link(png_filename,dpi=300):
    
    import matplotlib.pyplot as plt

    from IPython.display import display, FileLink, HTML

    # Save the plot to a file

    Save_png_name = "../PNG_Files/"+png_filename
    plt.savefig( Save_png_name , bbox_inches='tight', pad_inches=0.3, dpi=dpi )
    # Create a downloadable link
    #local_file = FileLink('./my_plot.png', result_html_prefix="Click here to download: ")
    #display(local_file)


    #local_file = FileLink('./my_plot.png', result_html_prefix="Click here to ")
    #display(HTML(f'{local_file._repr_html_()} download the plot'))

    # Assuming 'my_plot.png' is the file you want to download
    # Your custom filename

    # Create the HTML string, inserting the custom filename
    html_str = f'<a href="{Save_png_name}" download="{png_filename}">Download Plot as PNG</a>'
    # Display the download link
    display(HTML(html_str))