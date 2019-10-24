function [ cil ] = apply_offset_to_cilium( cil, offset_on_cilium, cilium_to_offset )
%apply_offset_to_cilium applies the offset calculated by
%find_offset_on_cilium to the correct cilium
%   Detailed explanation goes here


cnto = mod(cilium_to_offset,numel(cil)) + 1; % cilium not to offset

% apply offset on cto
cil(cilium_to_offset).tt = cil(cilium_to_offset).tt - cil(cilium_to_offset).tt(offset_on_cilium);
cil(cilium_to_offset).tt_um = cil(cilium_to_offset).tt_um - cil(cilium_to_offset).tt_um(offset_on_cilium); %now this has got 0 at offset


end

