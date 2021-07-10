function [ flux_out ] = route(flux_in, uh)

flux_out   = uh(1,1) * flux_in + uh(2,1);

end