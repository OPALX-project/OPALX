function [bz,gap_width]=get_magnet_bz(bz_data,max_gap_width, rad)

    bz = pchip(bz_data(1,:),bz_data(2,:), rad);
    gap_width = max_gap_width*bz_data(2,1)./bz;
end