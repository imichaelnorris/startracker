% recursive function to find root node for element 'tag' in equivalence 
% graph 'equiv' in a union-find fashion. Uses path compression for
% efficiency. equiv is a mapObj. Used as part of star segmentation in
% Star_Detector.m
function root = uf_root( tag, equiv )
    if equiv(tag) == tag
        root = tag;
    else
        equiv(tag) = uf_root(equiv(tag), equiv); % path compression
        root = equiv(tag);
    end
end