function R = limo_gae_weights(B, Brec, varargin)
%LIMO_GAE_WEIGHTS Score aligned GAE reconstructions without training or I/O.
% R = limo_gae_weights(B,Brec,'Name',value,...)
% Both inputs: beta x frame x channel x subject (singleton axes retained).
% Defaults: beta-space MAE, mean across betas, within-channel median/MAD,
% one-sided bisquare (threshold 3, width 4.685), raw weights in [0,1].
% See README.md for the complete option/output contract and limitations.
% No group labels, statistical tests, file writes, or implicit reordering.

if ~isnumeric(B) || ~isreal(B) || issparse(B) || isempty(B) || ndims(B)>4 || ...
        ~isnumeric(Brec) || ~isreal(Brec) || issparse(Brec) || ~isequal(size(B),size(Brec))
    error('limo_gae_weights:Shape','Inputs must be equal-sized real numeric arrays with at most four axes.');
end
d = [size(B,1),size(B,2),size(B,3),size(B,4)];
p = inputParser; p.FunctionName = mfilename; p.PartialMatching = false;
addParameter(p,'ErrorDomain','beta');
addParameter(p,'DesignMatrices',[]);
addParameter(p,'BetaIndices',1:d(1));
addParameter(p,'FrameIndices',1:d(2));
addParameter(p,'ChannelIndices',1:d(3));
addParameter(p,'Metric','mae');
addParameter(p,'BetaAggregation','mean');
addParameter(p,'Scope','channel');
addParameter(p,'Standardization','median-mad');
addParameter(p,'Reference','within-channel');
addParameter(p,'WeightFunction','bisquare');
addParameter(p,'Threshold',3);
addParameter(p,'TuningConstant',4.685);
addParameter(p,'Normalization','none');
addParameter(p,'MissingPolicy','error');
addParameter(p,'ReferenceMask',true(1,d(4)));
addParameter(p,'MinReferenceSubjects',3);
addParameter(p,'MinValidObservations',1);
addParameter(p,'ScaleTolerance',1e-12);
addParameter(p,'SubjectIDs',[]);
addParameter(p,'ChannelIDs',[]);
addParameter(p,'FrameValues',[]);
addParameter(p,'PredictionInfo',struct());
addParameter(p,'ReturnReconstruction',false);
parse(p,varargin{:}); o = p.Results;
names = {'ErrorDomain','Metric','BetaAggregation','Scope','Standardization', ...
    'Reference','WeightFunction','Normalization','MissingPolicy'};
choices = {{'beta','fitted'},{'mae','rmse'},{'mean','median','none'}, ...
    {'channel','subject','both'},{'median-mad','mean-sd','none'}, ...
    {'within-channel','pooled'},{'bisquare','powerlaw','none'}, ...
    {'none','mean-one'},{'error','omit'}};
for k=1:numel(names)
    name = names{k}; value = o.(name);
    if isstring(value) && isscalar(value), value = char(value); end
    if ~ischar(value) || size(value,1)~=1 || ~any(strcmpi(value,choices{k}))
        error('limo_gae_weights:Option','Invalid %s.',name);
    end
    o.(name) = lower(value);
end
o.BetaIndices = indices(o.BetaIndices,d(1));
o.FrameIndices = indices(o.FrameIndices,d(2));
o.ChannelIndices = indices(o.ChannelIndices,d(3));
scalarcheck(o.Threshold,0,Inf,false);
scalarcheck(o.TuningConstant,0,Inf,false);
if o.TuningConstant==0, error('limo_gae_weights:Option','TuningConstant must be positive.'); end
scalarcheck(o.ScaleTolerance,0,1,false);
scalarcheck(o.MinReferenceSubjects,2,Inf,true);
scalarcheck(o.MinValidObservations,1,Inf,true);
if ~strcmp(o.WeightFunction,'bisquare') && ~any(strcmp(p.UsingDefaults,'TuningConstant'))
    error('limo_gae_weights:Option','TuningConstant applies only to bisquare.');
end
if ~islogical(o.ReferenceMask) || ~isvector(o.ReferenceMask) || numel(o.ReferenceMask)~=d(4)
    error('limo_gae_weights:Option','ReferenceMask must be a logical vector, one entry per subject.');
end
o.ReferenceMask = reshape(o.ReferenceMask,1,[]);
if ~islogical(o.ReturnReconstruction) || ~isscalar(o.ReturnReconstruction) || ...
        ~isstruct(o.PredictionInfo) || ~isscalar(o.PredictionInfo)
    error('limo_gae_weights:Option','ReturnReconstruction must be logical scalar; PredictionInfo a scalar struct.');
end
checkids(o.SubjectIDs,d(4)); checkids(o.ChannelIDs,d(3));
if ~isempty(o.FrameValues) && (~isnumeric(o.FrameValues) || ~isreal(o.FrameValues) || ...
        ~isvector(o.FrameValues) || numel(o.FrameValues)~=d(2) || any(~isfinite(o.FrameValues)))
    error('limo_gae_weights:Metadata','FrameValues must be a finite numeric vector matching the input frames.');
end
if strcmp(o.ErrorDomain,'fitted')
    if ~strcmp(o.BetaAggregation,'mean')
        error('limo_gae_weights:Option','BetaAggregation is not applicable in fitted mode.');
    end
    if ~iscell(o.DesignMatrices) || numel(o.DesignMatrices)~=d(4)
        error('limo_gae_weights:Design','Provide one design matrix per subject.');
    end
    for s=1:d(4)
        X = o.DesignMatrices{s};
        if ~isnumeric(X) || ~isreal(X) || ~ismatrix(X) || isempty(X) || ...
                size(X,2)~=d(1) || any(~isfinite(X(:)))
            error('limo_gae_weights:Design','Each design must be finite, real, trial x original beta.');
        end
    end
elseif ~isempty(o.DesignMatrices)
    error('limo_gae_weights:Option','DesignMatrices applies only to fitted mode.');
end

B = double(B(o.BetaIndices,o.FrameIndices,o.ChannelIndices,:));
Brec = double(Brec(o.BetaIndices,o.FrameIndices,o.ChannelIndices,:));
bad = ~isfinite(B) | ~isfinite(Brec);
if strcmp(o.MissingPolicy,'error') && any(bad(:))
    error('limo_gae_weights:Missing','Selected data contain NaN or Inf.');
end
B(bad)=NaN; Brec(bad)=NaN;
nb=numel(o.BetaIndices); nf=numel(o.FrameIndices); nc=numel(o.ChannelIndices); ns=d(4);
ne=nb; if strcmp(o.ErrorDomain,'fitted'), ne=1; end
E=NaN(ne,nc,ns); counts=zeros(ne,nc,ns);
for s=1:ns
    for ch=1:nc
        if strcmp(o.ErrorDomain,'fitted')
            % A missing coefficient invalidates that frame's modeled response.
            delta=reshape(B(:,:,ch,s)-Brec(:,:,ch,s),nb,nf);
            delta=double(o.DesignMatrices{s}(:,o.BetaIndices))*delta;
            [E(1,ch,s),counts(1,ch,s)]=metric(delta,o);
        else
            for b=1:nb
                [E(b,ch,s),counts(b,ch,s)]=metric(B(b,:,ch,s)-Brec(b,:,ch,s),o);
            end
        end
    end
end
R.errors.by_beta_channel=[];
if strcmp(o.ErrorDomain,'beta'), R.errors.by_beta_channel=E; end
keep=strcmp(o.BetaAggregation,'none') && strcmp(o.ErrorDomain,'beta');
if keep
    A=E;
else
    A=NaN(1,nc,ns);
    for s=1:ns
        for ch=1:nc
            % Require every selected beta: do not change beta mix by subject.
            v=E(:,ch,s);
            if all(isfinite(v))
                if strcmp(o.BetaAggregation,'median'), A(1,ch,s)=median(v);
                else, A(1,ch,s)=mean(v); end
            end
        end
    end
end
na=size(A,1); Z=NaN(size(A)); ref=cell(na,nc);
for b=1:na
    if strcmp(o.Reference,'pooled')
        v=reshape(A(b,:,:),nc,ns);
        r=v(:,o.ReferenceMask);
        nref=sum(o.ReferenceMask & any(isfinite(v),1));
        [z,info]=standardize(v,r(:),nref,o);
        Z(b,:,:)=reshape(z,1,nc,ns);
        for ch=1:nc, ref{b,ch}=info; end
    else
        for ch=1:nc
            v=reshape(A(b,ch,:),1,ns); r=v(o.ReferenceMask);
            [z,ref{b,ch}]=standardize(v,r,sum(isfinite(r)),o);
            Z(b,ch,:)=reshape(z,1,1,ns);
        end
    end
end
R.errors.channel=channelshape(A,keep,nc,ns);
R.scores.channel=channelshape(Z,keep,nc,ns);
R.valid.channel=isfinite(R.scores.channel);
R.valid.observation_count=counts;
R.valid.selected_input=~bad;
R.reference.channel=ref;
R.errors.subject=NaN(na,ns);
S=NaN(na,ns);
for b=1:na
    for s=1:ns
        v=A(b,:,s); z=Z(b,:,s);
        % Require all selected channels for a comparable subject summary.
        if all(isfinite(v)), R.errors.subject(b,s)=mean(v); end
        if all(isfinite(z)), S(b,s)=mean(z); end
    end
end
R.weights=struct(); R.flags=struct(); R.diagnostics=struct();
if any(strcmp(o.Scope,{'channel','both'}))
    [q,w,f,ess]=weights(reshape(Z,na*nc,ns),o);
    R.weights.raw_channel=channelshape(reshape(q,na,nc,ns),keep,nc,ns);
    R.weights.channel=channelshape(reshape(w,na,nc,ns),keep,nc,ns);
    R.flags.channel.attenuated=channelshape(reshape(f,na,nc,ns),keep,nc,ns);
    R.flags.channel.zero_weight=R.weights.raw_channel==0;
    R.diagnostics.channel_effective_n=reshape(ess,na,nc);
end
if any(strcmp(o.Scope,{'subject','both'}))
    R.reference.subject=cell(na,1); SZ=NaN(na,ns);
    for b=1:na
        r=S(b,o.ReferenceMask);
        [SZ(b,:),R.reference.subject{b}]=standardize(S(b,:),r,sum(isfinite(r)),o);
    end
    R.scores.subject=SZ; R.valid.subject=isfinite(SZ);
    [q,w,f,ess]=weights(SZ,o);
    R.weights.raw_subject=q; R.weights.subject=w;
    R.flags.subject.attenuated=f; R.flags.subject.zero_weight=q==0;
    R.diagnostics.subject_effective_n=ess;
end
R.diagnostics.invalid_channel_scores=sum(~isfinite(Z(:)));
R.options=o;
R.provenance.prediction=o.PredictionInfo;
R.provenance.subject_ids=o.SubjectIDs;
R.provenance.channel_ids=o.ChannelIDs;
if ~isempty(o.ChannelIDs), R.provenance.channel_ids=o.ChannelIDs(o.ChannelIndices); end
R.provenance.frame_values=o.FrameValues;
if ~isempty(o.FrameValues), R.provenance.frame_values=o.FrameValues(o.FrameIndices); end
R.provenance.input_shape=d;
R.provenance.version='stage1-1';
R.provenance.channel_layout='channel x subject';
if keep, R.provenance.channel_layout='beta x channel x subject'; end
R.provenance.subject_layout='1 x subject';
if keep, R.provenance.subject_layout='beta x subject'; end
if o.ReturnReconstruction, R.reconstructed=Brec; end
end

function idx=indices(idx,n)
if ~isnumeric(idx) || ~isreal(idx) || ~isvector(idx) || isempty(idx) || ...
        any(~isfinite(idx)) || any(idx<1 | idx>n | idx~=fix(idx)) || numel(unique(idx))~=numel(idx)
    error('limo_gae_weights:Indices','Indices must be unique integers within the corresponding axis.');
end
idx=reshape(idx,1,[]);
end

function scalarcheck(x,lo,hi,integer)
if ~isnumeric(x) || ~isreal(x) || ~isscalar(x) || ~isfinite(x) || x<lo || x>hi || (integer && x~=fix(x))
    error('limo_gae_weights:Option','Invalid numeric option.');
end
end

function checkids(ids,n)
if isempty(ids), return; end
if ~(isnumeric(ids) || isstring(ids) || iscellstr(ids)) || ~isvector(ids) || numel(ids)~=n
    error('limo_gae_weights:Metadata','IDs must be a numeric/string/cellstr vector matching the input axis.');
end
if isnumeric(ids) && (~isreal(ids) || any(~isfinite(ids)))
    error('limo_gae_weights:Metadata','Numeric IDs must be finite and real.');
end
if numel(unique(ids))~=n
    error('limo_gae_weights:Metadata','IDs must be unique.');
end
end

function [e,n]=metric(delta,o)
if any(isinf(delta(:)))
    error('limo_gae_weights:Overflow','Finite input arithmetic overflowed; rescale the inputs.');
end
v=abs(delta(:)); v=v(isfinite(v)); n=numel(v); e=NaN;
if n<o.MinValidObservations, return; end
if strcmp(o.Metric,'mae'), e=mean(v);
else
    % Scale before squaring to avoid overflow for finite large errors.
    m=max(v); if m==0, e=0; else, e=m*sqrt(mean((v/m).^2)); end
end
end

function [z,info]=standardize(v,r,nref,o)
z=NaN(size(v)); r=r(isfinite(r)); good=isfinite(v);
info=struct('location',NaN,'scale',NaN,'status','insufficient-reference', ...
    'n_subjects',nref,'n_values',numel(r));
if strcmp(o.Standardization,'none')
    z(good)=v(good); info.location=0; info.scale=1; info.status='unstandardized'; return;
end
if nref<o.MinReferenceSubjects || isempty(r), return; end
if strcmp(o.Standardization,'mean-sd')
    loc=mean(r); scale=std(r,0); status='mean-sd';
else
    loc=median(r); scale=1.4826*median(abs(r-loc)); status='median-mad';
end
tol=o.ScaleTolerance*max(abs(r));
if scale<=tol && strcmp(o.Standardization,'median-mad')
    scale=(quantilelinear(r,0.75)-quantilelinear(r,0.25))/1.349;
    status='iqr-fallback';
end
info.location=loc; info.scale=scale;
if ~isfinite(loc) || ~isfinite(scale)
    info.status='nonfinite-scale'; return;
end
if scale<=tol
    if max(r)-min(r)<=tol && all(abs(v(good)-loc)<=tol)
        z(good)=0; info.status='constant';
    else
        info.status='degenerate-scale';
    end
else
    z(good)=(v(good)-loc)/scale; info.status=status;
end
end

function q=quantilelinear(x,p)
% Explicit interpolation at 1+(n-1)*p; no Statistics Toolbox needed.
x=sort(x(:)); h=1+(numel(x)-1)*p; a=floor(h); b=ceil(h);
q=x(a)+(h-a)*(x(b)-x(a));
end

function [q,w,flag,ess]=weights(z,o)
q=NaN(size(z)); good=isfinite(z); e=max(z(good)-o.Threshold,0);
switch o.WeightFunction
    case 'bisquare'
        v=zeros(size(e)); inside=e<o.TuningConstant;
        v(inside)=(1-(e(inside)/o.TuningConstant).^2).^2;
    case 'powerlaw'
        v=1./(1+e.^2);
    case 'none'
        v=ones(size(e));
end
q(good)=v; w=q; flag=good & z>o.Threshold;
if strcmp(o.WeightFunction,'none'), flag(:)=false; end
ess=NaN(size(z,1),1);
for row=1:size(q,1)
    valid=isfinite(q(row,:)); a=q(row,valid); den=sum(a);
    if isempty(a), continue; end
    if den==0
        error('limo_gae_weights:ZeroWeights','All valid raw weights in a row are zero.');
    end
    ess(row)=den^2/sum(a.^2);
    if strcmp(o.Normalization,'mean-one'), w(row,valid)=a*numel(a)/den; end
end
end

function v=channelshape(v,keep,nc,ns)
if ~keep, v=reshape(v,nc,ns); end
end
