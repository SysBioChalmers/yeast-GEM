function run_all()
% run_all  Rerun the analyses of the reviewer comments on one model and write
%   the tables and figures of the revision to results/.
%   run_all()          the model of the yeast-GEM checkout (YEAST_GEM_ROOT)
%
%   The flux sampling (r1_major5_degeneracy/runSampling.m) is run first as a
%   Slurm array with paper_figures/sample_array.sh, then run_all collects it.
%   Variants are named <release>_<constraints> (see common/variantModel.m); the
%   release is a yeast-GEM tag, or a variant added to loadTaggedModel for an
%   unreleased model.

setupStudy('3.0.0-rc1');
four = {'vdev_new','vdev_old','v902_new','v902_old'};

% R1 major 3: 2x2, undone components, opened constraints, curation replay
reproVariants(four, 'repro');
undone = {'undoTetraenol','undoMdh2','undoIdp2','undoAA','oldGAM','oldCofactor','mdh1old','fadh2x0'};
reproVariants([four, strcat('vdev_new+', undone)], 'function_ablation');
reproVariants({'vdev_new','vdev_new+gcvopen','vdev_new+glt1open','vdev_new+gcy1rev'}, 'sectionB');
curationAblation('curation_ablation', 'categories');
biomassTermVariants();

% R1 major 5 and R2 major 3: FVA on the optimal face and in the 1% band
for v = four
    fvaRun(v{1}, 0.99999);
    fvaRun(v{1}, 1.01);
end

% R2 major 1: deletion strains
knockoutPanel();
adh3Gpd2();
engineeredPanel();
gpdViability();
frdKnockouts();
for v = four, geneDeletionScreen(v{1}); end

% R2 major 3: proton budget
protonBudget();

% R1 major 5: tables and figures from the samples (after sample_array.sh)
if isfolder(fullfile(yeastRoot(), 'analysis', 'anaerobic', 'results', 'sampling'))
    analyzeSampling(0.99999);
    compareSensitivity();
    makeFigures(0.99999);
end
end
