% function to convert mouse to mouse number. 

function mouse_num = convertMouse(mouseID)

% check if the mouseID entered is a char type
%    if class(mouseID) == char
%        disp('converting...')
%    else
%        disp('check mouse ID')
%    end

% convert mouse 
% M433 = 1 
% M453 = 2
% M460 = 3
% M533 = 4 
% M534 = 5
% M545 = 6
% M547 = 7 
% M548 = 8
    
    if mouseID == ['M433']
        mouse_num = 1; 
    elseif mouseID == ['M453']
        mouse_num = 2; 
    elseif mouseID == ['M460']
        mouse_num = 3; 
    elseif mouseID == ['M533']
        mouse_num = 4; 
    elseif mouseID == ['M534']
        mouse_num = 5;
    elseif mouseID == ['M545']
        mouse_num = 6; 
    elseif mouseID == ['M547']
        mouse_num = 7; 
    else 
        mouse_num = 8; 
    end
    

end