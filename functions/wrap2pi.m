function [ wrapped_angle ] = wrap2pi( angle1 )
    %wrap2pi wraps radians from -pi to +pi around
    
    wrapped_angle = angle1;
    
    if (angle1 > pi) 
        wrapped_angle = -pi + (angle1-pi);
    end;
    
    if (angle1 < -pi) 
        wrapped_angle =  pi + (angle1+pi);
    end;
    
end

