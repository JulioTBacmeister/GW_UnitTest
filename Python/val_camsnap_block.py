###############################
# Validate camsnap 
###############################

tau_movmtn_gw=Xgw.TAU_MOVMTN.values
tau_movmtn_cam=X.TAU_MOVMTN.values
tau_rdg_gw=Xgw.TAU_RDG.values
tau_rdg_cam=X.TAU1RDGBETAM.values


fig,axs=plt.subplots( 1, 4 , figsize=(28,6) )

i=0
ax=axs[i]
ax.scatter( tau_movmtn_gw[0,20,:] , tau_movmtn_cam[0,20,:])
ax.set_title( 'moving mountain stress ' )

ax.set_xlabel( 'Unit Test'  )
ax.set_ylabel( 'CAM snapshot' )


i=1
ax=axs[1]
ax.scatter( tau_rdg_gw[0,20,:] , tau_rdg_cam[0,20,:])
ax.set_title( 'Ridge scheme stress ' )

ax.set_xlabel( 'Unit Test'  )
ax.set_ylabel( 'CAM snapshot' )


