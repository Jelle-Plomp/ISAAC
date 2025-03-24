function resetStartEvent()
% Function to reset the startEvent (this is different from jumping to the
% start, since that does not reset the startEvent)

% Change the startEvent back to the starting number for the live event
% If we don't do this, after unfreezing we would return to the event_nr_HFR
% meaning we save the data again.
    Control = evalin('base', 'Control');
    % This control_i prevents overwriting any present control arguments (if
    % this is not done, then some of the GUI buttons may fail)
    if isempty(Control.Command)
        control_i = 1;
    else
        control_i =2;
    end  
    event_nr_live = evalin('base', 'event_nr_live');
    Control(control_i).Command = 'set&Run';
    Control(control_i).Parameters = {'Parameters',1,'startEvent',event_nr_live};
    assignin('base','Control',Control);
end