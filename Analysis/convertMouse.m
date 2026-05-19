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
% M534 /= 5
% M452 = 5; 
% M545 = 6
% M547 = 7 
% M548 = 8
% M556 = 9
% M578 = 10 
% M600 = 11
% M646 = 12
% M648 = 13
% M650 = 14
% M654 = 15
    
    if mouseID == ['M433']
        mouse_num = 1; 
    elseif mouseID == ['M453']
        mouse_num = 2; 
    elseif mouseID == ['M460']
        mouse_num = 3; 
    elseif mouseID == ['M533']
        mouse_num = 4; 
    elseif mouseID == ['M452']
        mouse_num = 5;
    elseif mouseID == ['M545']
        mouse_num = 6; 
    elseif mouseID == ['M547']
        mouse_num = 7; 
    elseif mouseID == ['M548']
        mouse_num = 8; 
    elseif mouseID == ['M556']
        mouse_num = 9;
    elseif mouseID == ['M578']
        mouse_num = 10;
    elseif mouseID == ['M600']
        mouse_num = 11;
    elseif mouseID == ['M646']
        mouse_num = 12;
    elseif mouseID == ['M648']
        mouse_num = 13;
    elseif mouseID == ['M650']
        mouse_num = 14;
    elseif mouseID == ['M654']
        mouse_num = 15;
    end
    

end