function [fileNeedsToBeLoaded] = readDumps_valuesCheck(wvfld_filename, IT)
  % Test if the iteration field fills out its allowed range.
  itname_fills_number_span = regexp(wvfld_filename,['d',num2str(IT),'_']);
  if(isempty(itname_fills_number_span))
    itname_fills_number_span = 0;
  end
  
  % Test if the iteration number exists, with a zero right before (does not fill its allowed range).
  itname_has_zero_right_before = regexp(wvfld_filename,['0',num2str(IT),'_']);
  if(isempty(itname_has_zero_right_before))
    itname_has_zero_right_before = 0;
  end
  
  itname_has_something_else_right_before_and_before_text = regexp(wvfld_filename,['0+[1-9]0*',num2str(IT),'_']);
  if(isempty(itname_has_something_else_right_before_and_before_text))
    itname_has_something_else_right_before_and_before_text = 0;
  end
  
  fileNeedsToBeLoaded = ( itname_fills_number_span | (itname_has_zero_right_before & not(itname_has_something_else_right_before_and_before_text)) );
end

