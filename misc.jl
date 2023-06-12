

function str_to_num(s::String)
    subject_pattern = r"[^0-9.-]"
    eachmatch(subject_pattern, s)
end

"SWIP_009_S1.txt"[11]
