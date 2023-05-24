function [ v_LVLH ] = ECEI2LVLH( i, OM, om, theta, v_ECEI )

[  ~, v_pf ] = ECEI2pf_rotate( i, OM, om, v_ECEI' );

[  ~, v_LVLH ] = pf2LVLH_rotate( theta, v_pf );



end