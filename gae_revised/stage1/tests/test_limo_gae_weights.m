function test_limo_gae_weights()
% Basic numerical and contract tests. Run with stage1 on the MATLAB path.
tests={@identity,@metrics,@fitted,@aggregation,@mapping,@normalization, ...
    @permutations,@degenerate,@missing,@subject_scope,@validation,@no_side_effects};
for k=1:numel(tests)
    tests{k}(); fprintf('PASS %s\n',func2str(tests{k}));
end
fprintf('%d test groups passed.\n',numel(tests));
end

function identity()
B=reshape(1:72,2,3,3,4);
R=limo_gae_weights(B,B,'Scope','both','ReturnReconstruction',true);
near(R.errors.channel,zeros(3,4)); near(R.weights.channel,ones(3,4));
near(R.weights.subject,ones(1,4)); near(R.reconstructed,B);
assert(strcmp(R.reference.channel{1,1}.status,'constant'));
% Explicit singleton frame/channel dimensions must not be squeezed away.
R=limo_gae_weights(zeros(2,1,1,4),zeros(2,1,1,4),'BetaAggregation','none');
assert(isequal(size(R.weights.channel),[2 1 4]));
end

function metrics()
B=zeros(1,2,1,4); P=repmat(reshape([3 4],1,2,1,1),1,1,1,4);
R=limo_gae_weights(B,P,'Standardization','none','WeightFunction','none');
near(R.errors.channel,3.5*ones(1,4));
R=limo_gae_weights(B,P,'Metric','rmse','Standardization','none','WeightFunction','none');
near(R.errors.channel,sqrt(12.5)*ones(1,4));
assert(all(R.valid.observation_count(:)==2));
P=reshape([1 2 3 4],1,1,1,4);
R=limo_gae_weights(zeros(size(P)),P,'Standardization','mean-sd');
near(R.scores.channel,([-1.5 -.5 .5 1.5])/sqrt(5/3));
end

function fitted()
B=zeros(2,2,1,4); P=repmat([1 2;3 4],1,1,1,4);
X=repmat({[1 0;0 2;1 -1]},1,4);
R=limo_gae_weights(B,P,'ErrorDomain','fitted','DesignMatrices',X, ...
    'Standardization','none','WeightFunction','none');
% Response discrepancies: [1 2;6 8;-2 -2], abs sum 21, six observations.
near(R.errors.channel,3.5*ones(1,4)); assert(isempty(R.errors.by_beta_channel));
R=limo_gae_weights(B,P,'ErrorDomain','fitted','DesignMatrices',X, ...
    'BetaIndices',2,'Standardization','none','WeightFunction','none');
near(R.errors.channel,3.5*ones(1,4));
R=limo_gae_weights(B,P,'ErrorDomain','fitted','DesignMatrices',X, ...
    'BetaIndices',1,'Standardization','none','WeightFunction','none');
near(R.errors.channel,ones(1,4));
end

function aggregation()
B=zeros(3,1,1,4); P=repmat(reshape([1 2 9],3,1,1,1),1,1,1,4);
args={'Standardization','none','WeightFunction','none'};
R=limo_gae_weights(B,P,args{:}); near(R.errors.channel,4*ones(1,4));
R=limo_gae_weights(B,P,args{:},'BetaAggregation','median'); near(R.errors.channel,2*ones(1,4));
R=limo_gae_weights(B,P,args{:},'BetaAggregation','none');
near(R.errors.channel,repmat(reshape([1 2 9],3,1,1),1,1,4));
end

function mapping()
z=[0 3 5.3425 7.685 10]; B=zeros(1,1,1,5); P=reshape(z,1,1,1,5);
R=limo_gae_weights(B,P,'Standardization','none');
near(R.weights.raw_channel,[1 1 0.5625 0 0]);
assert(isequal(R.flags.channel.attenuated,logical([0 0 1 1 1])));
R=limo_gae_weights(B,P,'Standardization','none','WeightFunction','powerlaw');
near(R.weights.channel,1./(1+max(z-3,0).^2));
end

function normalization()
P=reshape([1 2 3 4 20],1,1,1,5); B=zeros(size(P));
A=limo_gae_weights(B,P); R=limo_gae_weights(B,P,'Normalization','mean-one');
near(A.weights.channel,[1 1 1 1 0]);
near(R.weights.channel,[1.25 1.25 1.25 1.25 0]);
near(R.weights.raw_channel,A.weights.channel);
assert(isequal(R.flags.channel.attenuated,A.flags.channel.attenuated));
near(R.diagnostics.channel_effective_n,4);
end

function permutations()
P=zeros(2,2,2,5);
for b=1:2
    for ch=1:2, P(b,:,ch,:)=repmat(reshape((b+ch)*[1 2 3 4 20],1,1,1,5),1,2,1,1); end
end
B=zeros(size(P)); A=limo_gae_weights(B,P,'Scope','both'); order=[5 2 4 1 3];
R=limo_gae_weights(B(:,:,[2 1],order),P(:,:,[2 1],order),'Scope','both');
near(R.weights.channel,A.weights.channel([2 1],order)); near(R.weights.subject,A.weights.subject(:,order));
R=limo_gae_weights(B,P,'FrameIndices',2,'ChannelIndices',2, ...
    'SubjectIDs',1:5,'ChannelIDs',{'A','B'},'FrameValues',[10 20]);
near(R.errors.channel,A.errors.channel(2,:)); assert(isequal(R.provenance.channel_ids,{'B'}));
near(R.provenance.frame_values,20);
% Within-channel standardization is invariant to positive channel scaling.
near(A.weights.channel(1,:),A.weights.channel(2,:));
% Pooled calibration has one reference shared by the channels.
P=reshape([1 10 2 20 3 30 4 40],1,1,2,4);
R=limo_gae_weights(zeros(size(P)),P,'Reference','pooled','Standardization','mean-sd');
v=[1 10 2 20 3 30 4 40];
near(R.scores.channel,(reshape(v,2,4)-mean(v))/std(v));
assert(isequal(R.reference.channel{1},R.reference.channel{2}));
end

function degenerate()
B=zeros(1,1,1,5);
R=limo_gae_weights(B,reshape([0 0 0 1 2],size(B)));
assert(strcmp(R.reference.channel{1}.status,'iqr-fallback'));
near(R.reference.channel{1}.scale,1/1.349);
R=limo_gae_weights(B,reshape([0 0 0 0 2],size(B)));
assert(strcmp(R.reference.channel{1}.status,'degenerate-scale'));
assert(all(isnan(R.weights.channel)));
R=limo_gae_weights(B,B,'ReferenceMask',logical([1 1 0 0 0]));
assert(strcmp(R.reference.channel{1}.status,'insufficient-reference'));
R=limo_gae_weights(B,reshape([0 0 0 1 2],size(B)), ...
    'ReferenceMask',logical([1 1 1 0 0]));
assert(strcmp(R.reference.channel{1}.status,'degenerate-scale'));
% No fixed epsilon: scaling the errors preserves nondegenerate z scores.
P=reshape([1 2 3 4 20],size(B)); A=limo_gae_weights(B,P);
R=limo_gae_weights(B,P*1e-15); near(R.scores.channel,A.scores.channel);
end

function missing()
B=zeros(1,2,1,4); P=repmat(reshape([1 2 3 4],1,1,1,4),1,2,1,1);
P(1,1,1,1)=NaN;
throws(@()limo_gae_weights(B,P),'limo_gae_weights:Missing');
R=limo_gae_weights(B,P,'MissingPolicy','omit'); near(R.errors.channel,[1 2 3 4]);
assert(R.valid.observation_count(1,1,1)==1);
R=limo_gae_weights(B,P,'MissingPolicy','omit','MinValidObservations',2);
assert(isnan(R.errors.channel(1))); assert(isnan(R.weights.channel(1)));
P(1,2,1,1)=Inf;
R=limo_gae_weights(B,P,'MissingPolicy','omit'); assert(isnan(R.weights.channel(1)));
end

function subject_scope()
P=reshape([1 2 3 4 20],1,1,1,5); B=zeros(size(P));
R=limo_gae_weights(B,P,'Scope','subject');
assert(~isfield(R.weights,'channel')); near(R.weights.subject,[1 1 1 1 0]);
% Unstandardized pooled scoring still supports a single subject explicitly.
R=limo_gae_weights(0,0,'Reference','pooled','Standardization','none');
near(R.weights.channel,1);
end

function validation()
B=zeros(1,2,1,4);
throws(@()limo_gae_weights(B,zeros(2,1,1,4)),'limo_gae_weights:Shape');
throws(@()limo_gae_weights(B,B,'BetaIndices',[1 1]),'limo_gae_weights:Indices');
throws(@()limo_gae_weights(B,B,'Threshold',-1),'limo_gae_weights:Option');
throws(@()limo_gae_weights(B,B,'WeightFunction','powerlaw','TuningConstant',2),'limo_gae_weights:Option');
throws(@()limo_gae_weights(B,B,'ErrorDomain','fitted'),'limo_gae_weights:Design');
throws(@()limo_gae_weights(B,B,'SubjectIDs',[1 1 2 3]),'limo_gae_weights:Metadata');
throws(@()limo_gae_weights(B,10+ B,'Standardization','none'),'limo_gae_weights:ZeroWeights');
throws(@()limo_gae_weights(realmax*ones(size(B)),-realmax*ones(size(B))), ...
    'limo_gae_weights:Overflow');
end

function no_side_effects()
before=pwd; state=rng; listing=stable_listing(); B=zeros(1,1,1,4);
limo_gae_weights(B,B);
assert(strcmp(pwd,before)); assert(isequal(rng,state)); assert(isequal(stable_listing(),listing));
end

function snapshot=stable_listing()
% Octave's extra statinfo includes access times, which can change on reads.
% Check contents/modification metadata, not access metadata or '.'/'..'.
files=dir;
files=files(~ismember({files.name},{'.','..'}));
snapshot={files.name; files.isdir; files.bytes; files.datenum};
end

function near(a,b)
assert(isequal(size(a),size(b)),'Unexpected output shape.');
assert(all(abs(a(:)-b(:))<=1e-10*max(1,max(abs(b(:))))),'Numerical mismatch.');
end

function throws(f,id)
try, f(); catch err, assert(strcmp(err.identifier,id),err.message); return; end
error('Expected error %s.',id);
end
