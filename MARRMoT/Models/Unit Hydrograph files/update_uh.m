function [ uh ] = update_uh(uh, flux_in)

uh(2,:)   = (uh(1,:) .* flux_in) + uh(2,:);
uh(2,:)   = circshift(uh(2,:),-1);
uh(2,end) = 0;

end