# Differential Workflow

The detailed differential workflow API page has moved to
[Differential Analysis Workflow](workflows/differential-analysis.md).

Use `DifferentialAnalysisWorkflow` from top-level `phospy`, and use
`DifferentialAnalysisRequest`, `DifferentialAnalysisConfig`,
`ExperimentalDesign`, `SampleDesignRecord`, and `Contrast` from `phospy.api`.

```python
from phospy import DifferentialAnalysisWorkflow
from phospy.api import DifferentialAnalysisRequest
```

`from phospy import DifferentialAnalysis` and
`from phospy.api import DifferentialAnalysis` are not supported public routes.
