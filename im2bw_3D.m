function [out] = im2bw_3D(out, level)
    dim3=size(out,3);
    for i=1:dim3
        out(:,:,i) = im2bw(out(:,:,i), level);
    end

end