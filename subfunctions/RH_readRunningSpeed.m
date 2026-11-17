function [vecSamples,vecRunningData,vecSyncData,sStream]=RH_readRunningSpeed(sStream,intRunningChanNI,intStimOnsetChanNI)
	%RH_readRunningSpeed Read real-time running speed (new SpikeGLX SDK). Syntax:
	%   [vecSamples,vecRunningData,vecSyncData,sStream]=RH_readRunningSpeed(sStream,intRunningChanNI,intStimOnsetChanNI)
	%
	%Port of Acquipix/Helpers/readRunningSpeed.m to the (js, ip) SpikeGLX API
	%(GetStreamSampleCount / Fetch(hSGL,js,ip,...)).
	%
	%inputs:
	% - sStream: structure with fields hSGL, intStreamNI (js; 0=nidq),
	%   intSubStreamNI (ip; 0 for nidq), intLastFetchNI, dblSampFreqNI, NI2V
	% - intRunningChanNI: NI channel connected to running wheel
	% - intStimOnsetChanNI: (optional) NI channel connected to screen diode
	%
	%outputs:
	% - vecSamples: NI sample numbers corresponding to vecRunningData/vecSyncData
	% - vecRunningData: raw running wheel data (empty if nothing could be fetched)
	% - vecSyncData: screen diode data (empty if channel was not supplied)
	% - sStream: updated structure
	%
	%Based on readRunningSpeed v1.1 by Jorrit Montijn; ported by Robin Haak

	%% get variables
	if nargin < 3 || isempty(intStimOnsetChanNI)
		intStimOnsetChanNI = [];
	end
	vecSamples = [];
	vecRunningData = [];
	vecSyncData = [];

	hSGL = sStream.hSGL;
	intJS = sStream.intStreamNI;
	if isfield(sStream,'intSubStreamNI'),intIP = sStream.intSubStreamNI;else,intIP = 0;end
	intLastFetchNI = sStream.intLastFetchNI;
	dblSampFreqNI = sStream.dblSampFreqNI;

	%stream buffer is ~8 s; keep max retrieval lower to be safe
	dblStreamBufferSizeSGL = 5;
	intMaxFetchNI = round(dblSampFreqNI*dblStreamBufferSizeSGL);
	intDownsampleNI = 1;

	%% determine fetch window
	intCurCountNI = GetStreamSampleCount(hSGL, intJS, intIP);
	if intLastFetchNI == 0 || (intLastFetchNI < (intCurCountNI - intMaxFetchNI))
		intRetrieveSamplesNI = round(1*dblSampFreqNI); %initial fetch: last 1 second
		intRetrieveSamplesNI = min(intCurCountNI-1,intRetrieveSamplesNI);
		intStartFetch = intCurCountNI - intRetrieveSamplesNI;
	else
		intStartFetch = intLastFetchNI - round(0.2*dblSampFreqNI);
		intStartFetch = max(intStartFetch,1);
		intRetrieveSamplesNI = intCurCountNI - intStartFetch;
	end
	if intRetrieveSamplesNI > intMaxFetchNI
		fprintf('<< WARNING >> Requested NI Fetch was %.1fs, max query is %.1fs\n',...
			intRetrieveSamplesNI/dblSampFreqNI,intMaxFetchNI/dblSampFreqNI);
		intRetrieveSamplesNI = intMaxFetchNI;
		intStartFetch = intCurCountNI - intRetrieveSamplesNI;
	end
	if intRetrieveSamplesNI <= 0
		return;
	end

	%% fetch
	vecChans = [intRunningChanNI intStimOnsetChanNI];
	try
		[matDataNI,intHeadCt] = Fetch(hSGL, intJS, intIP, intStartFetch, intRetrieveSamplesNI, vecChans, intDownsampleNI);
	catch ME
		warning([mfilename ':FetchFailed'],'NI Fetch failed [%s]: %s',ME.identifier,ME.message);
		return;
	end
	if isempty(matDataNI),return;end

	%SpikeGLX returns channels in ascending order, not in the requested order
	[~,vecReorder] = ismember(vecChans,sort(vecChans)); %column index of each requested channel

	%% update stream state
	sStream.intLastFetchNI = intCurCountNI;
	sStream.dblEphysTimeNI = intCurCountNI/dblSampFreqNI;

	%% collect outputs
	vecSamples = intHeadCt + (0:intDownsampleNI:(intDownsampleNI*(size(matDataNI,1)-1)));
	vecRunningData = matDataNI(:,vecReorder(1));
	if ~isempty(intStimOnsetChanNI)
		vecSyncData = sStream.NI2V*single(matDataNI(:,vecReorder(2)));
	end
end
