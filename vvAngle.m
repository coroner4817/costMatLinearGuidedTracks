function angle = vvAngle(vec1, vec2)
% vvAngle Returns the smallest angle between two vectors
% Note: Whilst relatively fast, the precision of this is limited near 0:
% vvAngle([1 1],[1 1]) is 1.2074e-06
    
    angle = abs(acosd(dot(vec1, vec2) / (norm(vec1) * norm(vec2))));
    
    % this should not be necessary, I think
    if angle > 180
        angle = 360 - angle;
    end
end

