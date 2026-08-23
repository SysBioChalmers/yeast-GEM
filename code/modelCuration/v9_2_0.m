% This scripts applies curations to be applied on yeast-GEM release 9.1.0,
% to get to yeast-GEM release 9.2.0.
%
% Curations in this release:
%   - Correction of erroneous metabolite annotations (Issue #378)
%   - Addition of a mitochondrial acetate transporter (Issue #388)
%   - Addition of ribonucleotide reductase acting on UDP (Issue #380)

%% Load yeast-GEM 9.1.0 (requires local yeast-GEM git repository)
cd ..
codeDir=pwd();
model = getEarlierModelVersion('9.1.0');
model.id='yeastGEM_develop';
model.version='';
% dataDir=fullfile(pwd(),'..','data','modelCuration','v9.2.0'); % No dataDir required for these curations
cd modelCuration

%% ========================================================================
% Correct erroneous metabolite annotations (Issue #378)
%
% Distinct metabolites were annotated with identical identifiers, and in
% several cases the identifier referred to an entirely different compound.
% Each identifier below was verified against the ChEBI, KEGG, BiGG and
% MetaNetX web services. Where a MetaNetX identifier had to be replaced,
% the current (non-deprecated) identifier is used; identifiers that are
% merely deprecated but still refer to the correct compound are left
% untouched, as their bulk update is covered by Issue #245.

% ITP[m] carried a complete set of succinate annotations, copied from
% succinate[m]. Align them with the already correct ITP[c] (s_0950), and
% add the MetaCyc and ModelSEED identifiers for ITP.
model = editMiriam(model,'met','s_4313','bigg.metabolite','itp','replace');
model = editMiriam(model,'met','s_4313','metanetx.chemical','MNXM423','replace');
model = editMiriam(model,'met','s_4313','biocyc','ITP','replace');            % was SUC = succinate
model = editMiriam(model,'met','s_4313','seed.compound','cpd00068','replace'); % was cpd00036 = succinate
model = removeMiriam(model,'met','s_4313','pubchem.compound');      % was CID 160419 = succinate

% L-xylulose[c] carried the ChEBI and MetaNetX identifiers of D-xylulose,
% which are correctly assigned to D-xylulose[c] (s_0580).
model = editMiriam(model,'met','s_4180','chebi','CHEBI:17399','replace');
model = editMiriam(model,'met','s_4180','metanetx.chemical','MNXM1371095','replace');

% The protein and carbohydrate biomass pseudometabolites carried
% identifiers of real compounds: bigg:protein is Torasemide-M3,
% kegg:C00492 and MNXM621 are raffinose (correctly assigned to raffinose,
% s_3715/s_4193), kegg:C05402 is melibiose and MNXM1434 is polydextrose.
model = removeMiriam(model,'met','s_3717','bigg.metabolite');
model = removeMiriam(model,'met','s_3717','kegg.compound');
model = removeMiriam(model,'met','s_3717','metanetx.chemical');
model = removeMiriam(model,'met','s_3718','kegg.compound');
model = removeMiriam(model,'met','s_3718','metanetx.chemical');
% The cofactor pseudometabolite carried the SMILES of calcium levofolinate
model.metSmiles(getIndexes(model,'s_4205','mets')) = {''};

% MNXM390 is L-galactose, while it was assigned to both D-galactose (in
% three compartments) and alpha-D-galactose. In addition, D-galactose[v]
% was annotated as D-galactopyranose (CHEBI:4139), unlike its counterparts
% in the cytoplasm and extracellular space.
model = editMiriam(model,'met',{'s_0558','s_0559','s_3780'},'metanetx.chemical','MNXM735266','replace');
model = editMiriam(model,'met','s_3780','chebi','CHEBI:12936','replace');
model = editMiriam(model,'met','s_3862','metanetx.chemical','MNXM1364486','replace');

% ribose-5-phosphate[c] has charge -2, while CHEBI:18189 is the neutral
% alpha-D-ribofuranose 5-phosphate. CHEBI:140679 is the (2-) form, which
% is also consistent with its bigg:r5p and kegg:C03736 annotations.
model = editMiriam(model,'met','s_1408','chebi','CHEBI:140679','replace');

% D-ribofuranose 5-phosphate[c] and aldehydo-D-ribose 5-phosphate[c] are
% distinct metabolites, but shared both their KEGG and MetaNetX
% identifiers. MetaNetX has a dedicated entry for the open-chain form,
% while kegg:C00117 is the unspecific D-ribose 5-phosphate that matches
% neither of the two in particular.
model = editMiriam(model,'met','s_3768','metanetx.chemical','MNXM1363910','replace');
model = removeMiriam(model,'met','s_3768','kegg.compound');
model = editMiriam(model,'met','s_3852','metanetx.chemical','MNXM1363911','replace');
model = removeMiriam(model,'met','s_3852','kegg.compound');

% 6-(N-acetyl-alpha-D-glucosaminyl)-1-phosphatidyl-1D-myo-inositol[er] and
% G00143[er] shared CHEBI:57265, which is the (1-) form and therefore only
% matches G00143[er] (charge -1). The former (charge 0) is CHEBI:12194.
model = editMiriam(model,'met','s_0331','chebi','CHEBI:12194','replace');
model = editMiriam(model,'met','s_0331','metanetx.chemical','MNXM1104777','replace');
model = editMiriam(model,'met','s_3879','metanetx.chemical','MNXM741477','replace');

% L-Methionine S-oxide is the unspecific sulfoxide, while L-methionine
% (S)-S-oxide is the (S)-epimer. Their BiGG identifiers were swapped
% (bigg:metsox_S__L is the (S)-epimer, bigg:metox is unspecific), and both
% MetaNetX identifiers resolve to the (S)-epimer.
model = editMiriam(model,'met','s_3807','bigg.metabolite','metsox_S__L','replace');
model = editMiriam(model,'met',{'s_3837','s_4176'},'bigg.metabolite','metox','replace');
model = editMiriam(model,'met',{'s_3837','s_4176'},'metanetx.chemical','MNXM1364200','replace');

% The names of taurocholate and taurocholic acid were swapped: the
% metabolites named "taurocholic acid" have formula C26H44NO7S and charge
% -1 and are annotated with CHEBI:36257 (taurocholate), while those named
% "taurocholate" have formula C26H45NO7S, charge 0 and CHEBI:28865
% (taurocholic acid). Rename to match formula, charge and ChEBI. Note that
% BiGG (tchola) and KEGG (C05122) do not distinguish between the two
% protonation states, so these identifiers legitimately remain shared.
model.metNames(getIndexes(model,{'s_1473','s_1474'},'mets')) = {'taurocholate'};
model.metNames(getIndexes(model,{'s_4267','s_4268'},'mets')) = {'taurocholic acid'};

%% ========================================================================
% Add a mitochondrial acetate transporter (Issue #388)
%
% Mitochondrial acetate is produced by Ach1 (r_0111) and the mitochondrial
% aldehyde dehydrogenases (r_0174, r_0175), and consumed by Acs1 (r_0113),
% but the model had no acetate transport across the mitochondrial
% membrane, so that cytosolic and mitochondrial acetate could only be
% connected via acetaldehyde diffusion. Ach1 transfers CoA from acetyl-CoA
% to succinate, after which the resulting acetate is exported to the
% cytosol (10.1093/femsyr/fov015). The transport mechanism is unknown, so
% no genes are assigned, matching the other unassigned acetate transport
% reactions (r_1635, r_4602).
% The reaction is irreversible in the export direction. Acetate import into
% the mitochondrion is reported for acetate detoxification
% (10.1016/j.fgb.2009.03.004), but that is a response to acetate stress,
% not a feature of growth on defined medium. Allowing import by default
% lets cytosolic acetate feed mitochondrial acetyl-CoA via Acs1 (r_0113),
% which shifts the anaerobic redox balance and degrades the predictions
% that release 9.1.0 was curated for: glycerol drops from 4.35 to 3.81
% (experimental 4.5 +/- 0.4), the mean relative error of the main
% fermentation products rises from 0.065 to 0.103, and the anaerobic flux
% prediction R2 drops from 0.9967 to 0.9950. Set the lower bound to -1000
% when explicitly simulating acetate stress.
newRxn = struct('rxns', {generateNewIds(model,'rxns','r_',1)}, ...
                'equations', {{'acetate[m] => acetate[c]'}}, ...
                'rxnNames', {{'acetate transport'}}, ...
                'subSystems', {{'Transport [c, m]'}}, ...
                'grRules', {{''}},...
                'rxnReferences', {{'10.1093/femsyr/fov015; 10.1016/j.fgb.2009.03.004'}}, ...
                'rxnMiriams', {{struct('name',{{'metanetx.reaction'}},'value',{{'MNXR95431'}})}}, ...
                'rxnNotes', {{'Set LB to -1000 to allow acetate import when simulating acetate stress'}}, ...
                'rxnConfidenceScores', {2});
model = addRxns(model,newRxn,3);
fprintf('Identifier of new reaction "%s": %s\n', newRxn.rxnNames{1}, model.rxns{end});

%% ========================================================================
% Add ribonucleotide reductase acting on UDP (Issue #380)
%
% The ribonucleotide reductase complex binds UDP (10.1073/pnas.0600443103)
% and reduces it to dUDP, as it does in other organisms
% (10.1016/S0076-6879(78)51032-X, 10.1242/dmm.050775). The reaction is
% present in the earlier S. cerevisiae models iMM904 and iND750 as RNDR4,
% but was missing from yeast-GEM, where dUDP could only be formed from
% dUMP. Added for both the cytoplasm and the nucleus, mirroring the
% existing reductases acting on ADP, CDP and GDP (r_0974-r_0979).
newRxns = struct('rxns', {generateNewIds(model,'rxns','r_',2)}, ...
                 'equations', {{'UDP[c] + H+[c] + TRX1[c] => dUDP[c] + H2O[c] + TRX1 disulphide[c]'; ...
                                'UDP[n] + H+[n] + TRX1[n] => dUDP[n] + H2O[n] + TRX1 disulphide[n]'}}, ...
                 'rxnNames', {{'ribonucleotide reductase';'ribonucleotide reductase'}}, ...
                 'subSystems', {{'Pyrimidine metabolism';'Pyrimidine metabolism'}}, ...
                 'grRules', {{'YER070W or YGR180C or YIL066C or YJL026W'; ...
                              'YER070W or YGR180C or YIL066C or YJL026W'}}, ...
                 'eccodes', {{'1.17.4.1';'1.17.4.1'}}, ...
                 'rxnReferences', {{'10.1073/pnas.0600443103';'10.1073/pnas.0600443103'}}, ...
                 'rxnConfidenceScores', [2;2]);
newRxns.rxnMiriams = {struct('name',{{'bigg.reaction';'kegg.reaction';'metanetx.reaction'; ...
                                      'kegg.pathway';'kegg.pathway';'kegg.pathway'}}, ...
                             'value',{{'RNDR4';'R02018';'MNXR104066'; ...
                                       'sce00230';'sce00240';'sce00480'}}); ...
                      struct('name',{{'bigg.reaction';'kegg.reaction';'metanetx.reaction'; ...
                                      'kegg.pathway';'kegg.pathway';'kegg.pathway'}}, ...
                             'value',{{'RNDR4n';'R02018';'MNXR104066'; ...
                                       'sce00230';'sce00240';'sce00480'}})};
model = addRxns(model,newRxns,3);
fprintf('Identifiers of new reactions "%s": %s and %s\n', newRxns.rxnNames{1}, model.rxns{end-1}, model.rxns{end});

%% ========================================================================

%% DO NOT CHANGE OR REMOVE THE CODE BELOW THIS LINE.
% Show some metrics:
cd(fullfile(codeDir,'modelTests'))
disp('Run gene essentiality analysis')
[new.accuracy,new.tp,new.tn,new.fn,new.fp] = essentialGenes(model);
fprintf('Genes in model: %d\n',numel(model.genes));
fprintf('Gene essentiality accuracy: %.4f\n', new.accuracy);
fprintf('True non-essential genes: %d\n', numel(new.tp));
fprintf('True essential genes: %d\n', numel(new.tn));
fprintf('False non-essential genes: %d\n', numel(new.fp));
fprintf('False essential genes: %d\n', numel(new.fn));
fprintf('\nRun growth analysis\n')
R2=growth(model);
fprintf('R2 of growth prediction: %.4f\n', R2);

% Save model:
cd(codeDir)
saveYeastModel(model)
cd modelCuration
