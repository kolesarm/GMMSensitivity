function WriteChecksum(obj, checksum_file, label)
%
%	Helper function writes estimation result to checksum file.
%
    diary(checksum_file);
    diary on
    disp(label);
    obj.Play();
    diary off
end