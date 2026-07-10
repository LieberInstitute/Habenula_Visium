initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new(
    "ReducedDimensionPlot",
    Type = "UMAP",
    XAxis = 1L,
    YAxis = 2L,
    FacetRowByColData = "cell_type",
    FacetColumnByColData = "cell_type",
    ColorByColumnData = "cell_type",
    ColorByFeatureNameAssay = "logcounts",
    ColorBySampleNameColor = "#FF0000",
    ShapeByColumnData = "cell_type",
    # SizeByColumnData = "age",
    TooltipColumnData = character(0),
    FacetRowBy = "None",
    FacetColumnBy = "None",
    ColorBy = "Column data",
    ColorByDefaultColor = "#000000",
    ColorByFeatureName = "MBP",
    ColorByFeatureSource = "---",
    ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "Br9090_1_MHb.2",
    ColorBySampleSource = "---",
    ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None",
    SizeBy = "None",
    SelectionAlpha = 0.1,
    ZoomData = numeric(0),
    BrushData = list(),
    VisualBoxOpen = FALSE,
    VisualChoices = "Color",
    ContourAdd = FALSE,
    ContourColor = "#0000FF",
    PointSize = 1,
    PointAlpha = 1,
    Downsample = FALSE,
    DownsampleResolution = 200,
    CustomLabels = FALSE,
    CustomLabelsText = "Br9090_1_MHb.2",
    FontSize = 1,
    LegendPointSize = 1,
    LegendPosition = "Bottom",
    HoverInfo = TRUE,
    LabelCenters = FALSE,
    LabelCentersBy = "cell_type",
    LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(
        list(c(2L, 10L, 0L)),
        class = c(
            "package_version",
            "numeric_version"
        )
    )),
    PanelId = c(ReducedDimensionPlot = 1L),
    PanelHeight = 600L,
    PanelWidth = 5L,
    SelectionBoxOpen = FALSE,
    RowSelectionSource = "---",
    ColumnSelectionSource = "---",
    DataBoxOpen = FALSE,
    RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE,
    SelectionHistory = list()
)

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new(
    "ComplexHeatmapPlot",
    Assay = "logcounts",
    CustomRows = TRUE,
    CustomRowsText = "MME\nGPR149\nMMRN1\nCRYM\nANKRD18B\nMYOZ3\nIPCEF1\nPLEKHD1\nTMEM52\nECEL1\nSLC22A6\nTHBD\nZNF710\nHOXD1\nCOL20A1",
    ClusterRows = FALSE,
    ClusterRowsDistance = "spearman",
    ClusterRowsMethod = "ward.D2",
    DataBoxOpen = FALSE,
    VisualChoices = "Annotations",
    ColumnData = "cell_type",
    RowData = character(0),
    CustomBounds = FALSE,
    LowerBound = NA_real_,
    UpperBound = NA_real_,
    AssayCenterRows = TRUE,
    AssayScaleRows = TRUE,
    DivergentColormap = "blue < white < orange",
    ShowDimNames = "Rows",
    LegendPosition = "Right",
    LegendDirection = "Vertical",
    VisualBoxOpen = FALSE,
    NamesRowFontSize = 10,
    NamesColumnFontSize = 10,
    ShowColumnSelection = TRUE,
    OrderColumnSelection = TRUE,
    VersionInfo = list(iSEE = structure(
        list(c(2L, 10L, 0L)),
        class = c(
            "package_version",
            "numeric_version"
        )
    )),
    PanelId = c(ComplexHeatmapPlot = 1L),
    PanelHeight = 600L,
    PanelWidth = 7L,
    SelectionBoxOpen = FALSE,
    RowSelectionSource = "---",
    ColumnSelectionSource = "---",
    RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE,
    SelectionHistory = list()
)

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new(
    "RowDataTable",
    Selected = "OPRM1",
    Search = "",
    SearchColumns = c(
        "",
        "", "", "", "", "", ""
    ),
    HiddenColumns = character(0),
    VersionInfo = list(iSEE = structure(
        list(c(2L, 10L, 0L)),
        class = c(
            "package_version",
            "numeric_version"
        )
    )),
    PanelId = c(RowDataTable = 1L),
    PanelHeight = 600L,
    PanelWidth = 5L,
    SelectionBoxOpen = FALSE,
    RowSelectionSource = "---",
    ColumnSelectionSource = "---",
    DataBoxOpen = FALSE,
    RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE,
    SelectionHistory = list()
)

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new(
    "FeatureAssayPlot",
    Assay = "logcounts",
    XAxis = "Column data",
    XAxisColumnData = "cell_type",
    XAxisFeatureName = "MBP",
    XAxisFeatureSource = "---",
    XAxisFeatureDynamicSource = FALSE,
    YAxisFeatureName = "OPRM1",
    YAxisFeatureSource = "RowDataTable1",
    YAxisFeatureDynamicSource = TRUE,
    FacetRowByColData = "cell_type",
    FacetColumnByColData = "cell_type",
    ColorByColumnData = "cell_type",
    ColorByFeatureNameAssay = "logcounts",
    ColorBySampleNameColor = "#FF0000",
    ShapeByColumnData = "cell_type",
    # SizeByColumnData = "age",
    TooltipColumnData = character(0),
    FacetRowBy = "None",
    FacetColumnBy = "None",
    ColorBy = "Column data",
    ColorByDefaultColor = "#000000",
    ColorByFeatureName = "MBP",
    ColorByFeatureSource = "---",
    ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "Br9090_1_MHb.2",
    ColorBySampleSource = "---",
    ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None",
    SizeBy = "None",
    SelectionAlpha = 0.1,
    ZoomData = numeric(0),
    BrushData = list(),
    VisualBoxOpen = FALSE,
    VisualChoices = "Color",
    ContourAdd = FALSE,
    ContourColor = "#0000FF",
    PointSize = 1,
    PointAlpha = 1,
    Downsample = FALSE,
    DownsampleResolution = 200,
    CustomLabels = FALSE,
    CustomLabelsText = "Br9090_1_MHb.2",
    FontSize = 1,
    LegendPointSize = 1,
    LegendPosition = "Bottom",
    HoverInfo = TRUE,
    LabelCenters = FALSE,
    LabelCentersBy = "cell_type",
    LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(
        list(c(2L, 10L, 0L)),
        class = c("package_version", "numeric_version")
    )),
    PanelId = c(FeatureAssayPlot = 1L),
    PanelHeight = 600L,
    PanelWidth = 7L,
    SelectionBoxOpen = FALSE,
    RowSelectionSource = "---",
    ColumnSelectionSource = "---",
    DataBoxOpen = FALSE,
    RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE,
    SelectionHistory = list()
)
