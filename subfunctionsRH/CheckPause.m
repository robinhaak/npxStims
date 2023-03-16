function boolPause = CheckPause()
	%CheckEnter Checks if the Space + 'P' keys are pressed
	%   boolEnter=any(strcmpi(KbName(keyCode),'escape'));
	KbName('UnifyKeyNames');
	[keyIsDown, secs, keyCode] = KbCheck();
	boolPause=any(strcmpi(KbName(keyCode),'p')) & any(strcmpi(KbName(keyCode),'space'));              
end
