function plan = buildfile
plan = buildplan(localfunctions);

pkg_root = fullfile(plan.RootFolder, '+navsu');
test_root = fullfile(plan.RootFolder, 'navsu-unit-testing');
reportDir = fullfile(plan.RootFolder, 'reports');

%% self-tests
plan('test') = matlab.buildtool.tasks.TestTask(test_root, SourceFiles=pkg_root, IncludeSubfolders=true);

%% Coverage
if ~isMATLABReleaseOlderThan('R2024a')

  coverageReport = fullfile(reportDir, 'coverage-report.html');
  try
    report = matlabtest.plugins.codecoverage.StandaloneReport(coverageReport);
  catch
    report = coverageReport;
  end
  plan('coverage') = plan('test').addCodeCoverage(report);
  plan('coverage').DisableIncremental = true;
end
%% clean
if ~isMATLABReleaseOlderThan("R2023b")
  plan("clean") = matlab.buildtool.tasks.CleanTask();
end

end


function checkTask(context)
root = context.Plan.RootFolder;

c = codeIssues(root, IncludeSubfolders=true);

if isempty(c.Issues)
  fprintf('%d files checked OK with %s under %s\n', numel(c.Files), c.Release, root)
else
  disp(c.Issues)
  error("Errors found in " + join(c.Issues.Location, newline))
end

end
