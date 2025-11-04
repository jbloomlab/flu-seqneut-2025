import marimo

__generated_with = "0.17.6"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # Make summary plots of the titers of different strains against human sera
    """
    )
    return


@app.cell
def _():
    # Load context from pickled file.
    #
    # This cell supports multiple ways to provide context:
    # 1. Via command-line: marimo export html notebook.py -- --context-pickle path/to/context.pickle
    # 2. Via saved pickle: Manually save a context pickle to a dev location
    # 3. Stub context: If no pickle available, creates minimal empty context for exploration
    #
    # For interactive development with `marimo edit`, you can:
    # - Run the pipeline once to generate a real context pickle, then copy it to a dev location
    # - Or work with the stub context (downstream cells will show warnings/empty data)

    import argparse
    import os
    import pathlib
    import pickle
    import sys

    import marimo as mo

    # Check if context-pickle argument is provided (run by driver script)
    from_cmdline = "--context-pickle" in sys.argv

    if from_cmdline:
        # Running via driver script - parse args
        print("Loading context from command-line argument")
        p = argparse.ArgumentParser()
        p.add_argument("--context-pickle", required=True)
        args = p.parse_args()
        context_pickle_path = pathlib.Path(args.context_pickle)
    else:
        # Running in marimo edit - try to use development pickle
        print("Running in marimo edit mode")
        # set `context_pickle_path` to valid pickle if running via marimo edit
        context_pickle_path = None
        context_pickle_path = pathlib.Path(
            "results/aggregated_analyses/plot_human_sera_titers_context.pickle"
        )

    # Load context if pickle path exists and is valid
    if context_pickle_path and context_pickle_path.exists():
        print(f"Reading context from {context_pickle_path}")
        with open(context_pickle_path, "rb") as f_context:
            context = pickle.load(f_context)

        # Handle working directory
        context_workdir = context["workdir"]
        current_workdir = os.getcwd()

        if from_cmdline:
            # Running via snakemake - verify workdir matches
            if context_workdir != current_workdir:
                raise RuntimeError(
                    f"Context workdir mismatch!\n"
                    f"  Context was created in: {context_workdir}\n"
                    f"  Currently running in:   {current_workdir}\n"
                    f"This should not happen when running via Snakemake."
                )
            print(f"Verified working directory: {current_workdir}")
        else:
            # Running in marimo edit - change to context workdir
            if context_workdir and context_workdir != current_workdir:
                print(f"Changing directory from {current_workdir} to {context_workdir}")
                os.chdir(context_workdir)
            elif context_workdir:
                print(f"Already in correct working directory: {context_workdir}")
    else:
        # Create a minimal stub context for interactive development
        print("Creating minimal stub context that you need to complete")
        context = {
            "input": {},
            "output": {},
            "params": {},
            "wildcards": {},
            "threads": 1,
            "resources": {},
        }
    return context, mo


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Read data
    """
    )
    return


@app.cell
def _():
    import itertools
    import json

    import altair as alt

    import pandas as pd

    _ = alt.data_transformers.disable_max_rows()
    return alt, itertools, json, pd


@app.cell
def _(context):
    # Extract variables from context - raises KeyError if required keys missing
    metadata_csv = context["input"]["metadata_csv"]
    titers_csv = context["input"]["titers_csv"]
    virus_csv = context["input"]["virus_csv"]
    viral_strain_plot_order_csv = context["input"]["viral_strain_plot_order"]
    recent_vaccine_strains = context["params"]["recent_vaccine_strains"]
    circulating_strain_type = context["params"]["circulating_strain_type"]
    human_sera_plots_params = context["params"]["human_sera_plots_params"]
    summarized_titers_csv = context["output"]["summarized_titers_csv"]
    sera_collection_date_and_age_plot = context["output"][
        "sera_collection_date_and_age_plot"
    ]
    chart_htmls = context["output"]["chart_htmls"]
    return (
        chart_htmls,
        circulating_strain_type,
        human_sera_plots_params,
        metadata_csv,
        recent_vaccine_strains,
        sera_collection_date_and_age_plot,
        summarized_titers_csv,
        titers_csv,
        viral_strain_plot_order_csv,
        virus_csv,
    )


@app.cell
def _(
    circulating_strain_type,
    metadata_csv,
    mo,
    pd,
    recent_vaccine_strains,
    titers_csv,
    viral_strain_plot_order_csv,
    virus_csv,
):
    metadata_all = pd.read_csv(metadata_csv)
    mo.output.append(mo.md(f"Read {len(metadata_all)=} for sera from {metadata_csv=}"))

    titers_all = pd.read_csv(titers_csv)
    mo.output.append(mo.md(f"\nRead {len(titers_all)=} titers from {titers_csv=}"))
    assert set(titers_all["serum"]) == set(metadata_all["serum"])

    viruses_all = (
        pd.read_csv(virus_csv)[
            ["strain", "subtype", "strain_type", "subclade", "vaccine_type"]
        ]
        .drop_duplicates()
        .rename(columns={"strain": "virus"})
    )
    assert set(titers_all["virus"]).issubset(viruses_all["virus"])
    assert set(recent_vaccine_strains).issubset(viruses_all["virus"])
    if not set(viruses_all["strain_type"]).issubset(
        {circulating_strain_type, "vaccine"}
    ):
        raise ValueError(
            f"{viruses_all['strain_type']=} != ['vaccine', {circulating_strain_type=}]"
        )
    viruses_all["strain_type"] = viruses_all["strain_type"].where(
        ~viruses_all["virus"].isin(recent_vaccine_strains), "recent_vaccine"
    )
    mo.output.append(mo.md(f"\nRead {len(viruses_all)=} viruses from {virus_csv=}"))

    viral_strain_plot_order = pd.read_csv(viral_strain_plot_order_csv)[
        "strain"
    ].tolist()
    assert set(viruses_all["virus"]).issubset(viral_strain_plot_order)
    return metadata_all, titers_all, viral_strain_plot_order, viruses_all


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Get fraction of sera and viruses with titers
    Look at the fraction of sera and viruses with titers.
    We generally may want to drop sera that lack titers for many viruses, and viruses that lack titers for many sera.
    Depending on `min_frac_action` (specified in configuration), we either raise an error if any sera or viruses are below the fractions specified in the configuration, or drop any sera or titers below these fractions.
    In general for production runs you may want to raise an error and filter these at an upstream step.
    """
    )
    return


@app.cell
def _(
    alt,
    human_sera_plots_params,
    metadata_all,
    mo,
    pd,
    titers_all,
    viruses_all,
):
    min_frac_action = human_sera_plots_params["min_frac_action"]

    titers = titers_all.copy()
    metadata = metadata_all.copy()
    viruses = viruses_all.copy()

    for min_frac_type, vals, frac_vals, col, frac_var in [
        (
            "min_frac_strains",
            set(metadata["serum"]),
            set(viruses["virus"]),
            "serum",
            "virus",
        ),
        (
            "min_frac_sera",
            set(viruses["virus"]),
            set(metadata["serum"]),
            "virus",
            "serum",
        ),
    ]:
        min_frac = human_sera_plots_params[min_frac_type]
        mo.output.append(
            mo.md(f"\nChecking for {min_frac_type=} at cutoff of {min_frac=}")
        )
        frac_df = (
            titers.groupby(col, as_index=False)
            .aggregate(n_w_titers=pd.NamedAgg(col, "count"))
            .assign(frac_w_titers=lambda x: x["n_w_titers"] / len(frac_vals))
        )
        frac_df = pd.concat(
            [
                frac_df,
                pd.DataFrame(
                    {
                        col: [v for v in vals if v not in set(frac_df[col])],
                        "n_w_titers": 0,
                        "frac_w_titers": 0.0,
                    }
                ),
            ],
            ignore_index=True,
        )
        assert len(frac_df) == len(vals)
        assert (frac_df["frac_w_titers"] <= 1).all()

        frac_df["below_cutoff"] = frac_df["frac_w_titers"] < min_frac

        frac_chart = (
            alt.Chart(frac_df)
            .encode(
                alt.X(
                    "frac_w_titers",
                    title=f"fraction {frac_var} with titers",
                    scale=alt.Scale(domain=[0, 1]),
                ),
                alt.Y(col, sort=alt.SortField("frac_w_titers", order="descending")),
                alt.Fill(
                    "below_cutoff",
                    title=f"below cutoff of {min_frac_type} = {min_frac}",
                ),
                tooltip=[col, "n_w_titers", alt.Tooltip("frac_w_titers", format=".3f")],
            )
            .mark_bar()
            .properties(
                height=alt.Step(11),
                width=135,
                title=f"{frac_var} with titers for each {col}",
            )
            .configure_axis(labelLimit=500)
            .configure_legend(titleLimit=500)
        )

        mo.output.append(frac_chart)

        failed = (
            frac_df.query("below_cutoff")
            .sort_values("frac_w_titers")
            .reset_index(drop=True)
        )
        below_frac = set(failed[col])
        mo.output.append(
            mo.md(
                f"Overall, {len(below_frac)} {col} have titers for less than {min_frac} {frac_var}"
            )
        )
        if len(failed):
            mo.output.append(failed)

        if min_frac_action == "raise":
            if not failed.empty:
                raise ValueError(
                    frac_df.query("below_cutoff")
                    .sort_values("frac_w_titers")
                    .reset_index(drop=True)
                )
        elif min_frac_action == "drop":
            if not failed.empty:
                mo.output.append(mo.md(f"Dropping these {col}"))
                titers = titers[~titers[col].isin(below_frac)]
                if col == "serum":
                    metadata = metadata[~metadata[col].isin(below_frac)]
                if col == "virus":
                    viruses = viruses[~viruses[col].isin(below_frac)]
        else:
            raise ValueError(f"invalid {min_frac_action=}")

    mo.output.append(
        mo.md(
            f"\nAfter any filtering:"
            f"\n {len(titers)=} / {len(titers_all)=}"
            f"\n {len(viruses)=} / {len(viruses_all)=}"
            f"\n {len(metadata)=} / {len(metadata_all)=}"
        )
    )
    return metadata, titers, viruses


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Plot age and collection distributions of sera
    """
    )
    return


@app.cell
def _(alt, metadata, mo, sera_collection_date_and_age_plot):
    sera_base = (
        alt.Chart(
            metadata[["group", "serum", "collection_date_numerical", "age_years"]]
            .assign(
                n=lambda x: x.groupby("group")["serum"].transform("count"),
                group=lambda x: x["group"] + " (n=" + x["n"].astype(str) + ")",
            )
            .drop(columns="n")
        )
        .encode(alt.Row("group", title=None))
        .mark_bar()
    )

    sera_date_chart = (
        sera_base.encode(
            alt.X(
                "yearmonth(collection_date_numerical)",
                title="collection date",
                axis=alt.Axis(format="%b-%Y", labelAngle=270),
                scale=alt.Scale(nice=False, padding=3),
            ),
            alt.Y(
                "count()",
                title="number of sera",
                scale=alt.Scale(nice=False, padding=3),
            ),
        )
        .resolve_scale(y="independent")
        .properties(height=80, width=70)
    )

    sera_age_chart = (
        sera_base.encode(
            alt.X(
                "age_years",
                title="age (years)",
                bin=alt.Bin(step=5, anchor=0),
                scale=alt.Scale(nice=False, padding=3),
            ),
            alt.Y(
                "count()",
                title="number of sera",
                scale=alt.Scale(nice=False, padding=3),
            ),
        )
        .resolve_scale(y="independent")
        .properties(height=80, width=150)
    )

    sera_chart = (
        alt.hconcat(sera_date_chart, sera_age_chart)
        .configure_axis(grid=False, titleFontWeight="normal")
        .configure_header(
            title=None, labelOrient="top", labelFontSize=11, labelPadding=2
        )
        .configure_facet(spacing=7)
        .configure_view(stroke="black")
        .properties(
            title=alt.TitleParams(
                "collection dates and subject ages for sera",
                anchor="middle",
            ),
        )
    )

    mo.output.append(sera_chart)

    mo.output.append(mo.md(f"Saving chart to {sera_collection_date_and_age_plot}"))
    sera_chart.save(sera_collection_date_and_age_plot)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Plot all the titers

    ### Assign label colors by subclade / vaccine type
    Define color mapping from subclade (circulating strains) or vaccine type (vaccine strains) to colors for label coloring, then create an expression that can be passed to `altair` *labelColor*:
    """
    )
    return


@app.cell
def _(
    circulating_strain_type,
    human_sera_plots_params,
    json,
    mo,
    pd,
    viral_strain_plot_order,
    viruses,
):
    strain_color_prop = viruses.assign(
        strain=lambda x: pd.Categorical(
            x["virus"], viral_strain_plot_order, ordered=True
        ),
        color_prop=lambda x: x["subclade"].where(
            x["strain_type"] == circulating_strain_type, x["vaccine_type"] + " vaccine"
        ),
    ).sort_values("strain")

    assert strain_color_prop["color_prop"].notnull().all()
    assert set(viruses["virus"]) == set(strain_color_prop["strain"])

    viruses["color_prop"] = viruses["virus"].map(
        strain_color_prop.set_index("strain")["color_prop"].to_dict()
    )

    prop_colors = human_sera_plots_params["prop_colors"]

    other_prop_colors = human_sera_plots_params["other_prop_colors"]

    for _subtype in strain_color_prop["subtype"].unique():
        subtype_color_props = (
            strain_color_prop[strain_color_prop["subtype"] == _subtype]["color_prop"]
            .unique()
            .tolist()
        )
        props_not_yet_colored = [p for p in subtype_color_props if p not in prop_colors]
        if len(props_not_yet_colored) > len(other_prop_colors):
            raise ValueError(
                f"props_not_yet_colored={props_not_yet_colored!r} longer than other_prop_colors={other_prop_colors!r}"
            )
        prop_colors.update(dict(zip(props_not_yet_colored, other_prop_colors)))

    mo.output.append(
        pd.Series(prop_colors).rename("color").rename_axis("property").to_frame()
    )
    assert set(strain_color_prop["color_prop"]).issubset(prop_colors)

    strain_color_prop = strain_color_prop.assign(
        color=lambda x: x["color_prop"].map(prop_colors)
    )

    color_mapping = strain_color_prop.set_index("virus")["color"].to_dict()

    # make a different color map for each subtype as they are plotted separately
    labelColor_expr = f"({json.dumps(color_mapping)})[datum.label] || 'black'"
    return labelColor_expr, strain_color_prop


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Now make nicely formatted charts
    """
    )
    return


@app.cell
def _(
    alt,
    human_sera_plots_params,
    labelColor_expr,
    metadata,
    mo,
    pd,
    titers,
    viral_strain_plot_order,
    viruses,
):
    # First set up the base chart and selections

    facet_size = human_sera_plots_params["facet_size"]
    facet_orientation = human_sera_plots_params["facet_orientation"]
    assert facet_orientation in {"vertical", "horizontal"}, f"{facet_orientation=}"
    titer_encoding = "x" if facet_orientation == "vertical" else "y"

    assert len(titers) == len(titers.groupby(["serum", "virus"]))
    assert viruses["virus"].nunique() == len(
        viruses.groupby(["virus", "subtype", "strain_type"])
    )

    # for group labels within altair we calculate these to have n too
    groups = (
        pd.concat([metadata, metadata.assign(group="All")])
        .groupby("group", as_index=False)
        .aggregate(n=pd.NamedAgg("serum", "nunique"))
        .assign(group=lambda x: x["group"] + " (n=" + x["n"].astype(str) + ")")["group"]
        .tolist()
    )
    groups = ["All", *sorted(metadata["group"].unique())]
    mo.output.append(mo.md(f"Plotting {groups=}"))

    virus_selection = alt.selection_point(
        fields=["virus"], on="mouseover", empty=False, clear="mouseout", nearest=False
    )

    serum_selection = alt.selection_point(
        fields=["serum"], on="mouseover", empty=False, clear="mouseout", nearest=False
    )

    group_selection = alt.selection_point(
        fields=["group"], bind="legend", empty="all", toggle="true", clear=False
    )

    # select by color used to color strain labels
    color_prop_selection = alt.selection_point(
        fields=["color_prop"], bind="legend", empty="all", toggle="true", clear=False
    )

    max_age = 5 * int(metadata["age_years"].max() // 5) + 5
    assert all(metadata["age_years"] <= max_age)
    min_age_slider = alt.param(
        value=0,
        bind=alt.binding_range(
            min=0, max=max_age, step=5, name="minimum subject age (years)"
        ),
    )
    max_age_slider = alt.param(
        value=max_age,
        bind=alt.binding_range(
            min=0, max=max_age, step=5, name="maximum subject age (years)"
        ),
    )

    # make the chart base, using transform_lookup to make it as small as possible
    # by looking up serum-specific and virus-specific annotations
    titers_base_nolookup = (
        alt.Chart(titers[["serum", "virus", "titer"]])
        .add_params(
            virus_selection,
            serum_selection,
            group_selection,
            color_prop_selection,
            min_age_slider,
            max_age_slider,
        )
        .encode(
            **{
                ("y" if facet_orientation == "vertical" else "x"): alt.Y(
                    "virus",
                    sort=list(reversed(viral_strain_plot_order)),
                    axis=alt.Axis(
                        labelLimit=500,
                        labelColor={"expr": labelColor_expr},
                        labelFontWeight=600,  # make a bit bolder so colors show
                        labelExpr="replace(datum.label, regexp('_[^_]*$'), '')",  # remove _H1N1 or _H3N2
                    ),
                ),
            },
        )
        .properties(
            **(
                {"height": alt.Step(11), "width": facet_size}
                if facet_orientation == "vertical"
                else {"width": alt.Step(11), "height": facet_size}
            )
        )
    )

    # dummy chart to bind the selectable legend for serum group
    dummy_group_chart = (
        alt.Chart(pd.DataFrame({"group": groups}))
        .add_params(group_selection)
        .mark_point(opacity=0)
        .encode(
            fill=alt.Fill(
                "group",
                title="serum group (click to select)",
                scale=alt.Scale(domain=groups, range=["gray"]),
                legend=alt.Legend(
                    symbolStrokeColor="black", symbolOpacity=1, columns=6
                ),
            )
        )
        .properties(width=1, height=1)  # tiny plot; legend renders outside
    )

    # because of scoping issues when layering and faceting charts with
    # transform_lookups (faceting must be done before lookups), we add
    # this function to do the faceting and lookups
    def facet_and_add_lookups(chart):
        return (
            chart
            # facet
            .facet(
                {
                    (
                        "column" if facet_orientation == "vertical" else "row"
                    ): alt.Column("group_n:N", title=None)
                }
            )
            # lookup additional data
            .transform_lookup(
                lookup="serum",
                from_=alt.LookupData(
                    data=metadata,
                    key="serum",
                    fields=[
                        "group",
                        "collection_date_string",
                        "age_string",
                        "age_years",
                        "sex",
                    ],
                ),
            )
            .transform_lookup(
                lookup="virus",
                from_=alt.LookupData(
                    data=viruses,
                    key="virus",
                    fields=["subtype", "strain_type", "subclade", "color_prop"],
                ),
            )
            # trick to make a new variable with all groups
            .transform_calculate(facet_with_all="[datum.group, 'All']")
            .transform_flatten(["facet_with_all"], as_=["group"])
            # filter by property used to color strain labels
            .transform_filter(color_prop_selection)
            # filter by group and age
            .transform_filter(group_selection)
            .transform_filter(alt.datum["age_years"] >= min_age_slider)
            .transform_filter(alt.datum["age_years"] <= max_age_slider)
            # make facet labels w n per group
            .transform_joinaggregate(n_per_group="distinct(serum)", groupby=["group"])
            .transform_calculate(
                group_n="datum.group + ' (n=' + datum.n_per_group + ')'"
            )
        )

    return (
        color_prop_selection,
        dummy_group_chart,
        facet_and_add_lookups,
        facet_orientation,
        serum_selection,
        titer_encoding,
        titers_base_nolookup,
        virus_selection,
    )


@app.cell
def _(alt, human_sera_plots_params, mo):
    # set titer scale
    titer_lower_limit = human_sera_plots_params["titer_lower_limit"]
    mo.output.append(mo.md(f"Using {titer_lower_limit=}"))
    titer_scale = alt.Scale(
        type="log", nice=False, domainMin=titer_lower_limit, padding=4
    )
    return titer_lower_limit, titer_scale


@app.cell
def _(alt, titer_encoding, titer_scale, titers_base_nolookup, virus_selection):
    # make median titer point chart
    median_points = (
        titers_base_nolookup.transform_aggregate(
            median_titer="median(titer)",
            groupby=["virus", "subtype", "strain_type", "subclade", "group"],
        )
        .encode(
            **{
                titer_encoding: alt.X(
                    "median_titer:Q", title="titer", scale=titer_scale
                )
            },
            tooltip=[
                "virus",
                alt.Tooltip("median_titer:Q", format=".1f"),
                "strain_type:N",
                "subclade:N",
            ],
            color=alt.condition(virus_selection, alt.value("red"), alt.value("black")),
            size=alt.condition(virus_selection, alt.value(80), alt.value(40)),
        )
        .mark_circle(opacity=1)
    )

    # facet_and_add_lookups(median_points)
    return (median_points,)


@app.cell
def _(alt, serum_selection, titer_encoding, titer_scale, titers_base_nolookup):
    # make per-serum lines
    serum_lines = titers_base_nolookup.encode(
        **{titer_encoding: alt.X("titer", scale=titer_scale)},
        detail=alt.Detail("serum"),
        tooltip=[
            "virus",
            "serum",
            alt.Tooltip("titer", format=".1f"),
            alt.Tooltip("collection_date_string:N", title="serum date"),
            alt.Tooltip("age_string:N", title="age"),
            "sex:N",
        ],
        size=alt.condition(serum_selection, alt.value(3), alt.value(1.5)),
        opacity=alt.condition(serum_selection, alt.value(1), alt.value(0.2)),
    ).mark_line()

    # facet_and_add_lookups(serum_lines + median_points)
    return (serum_lines,)


@app.cell
def _(alt, titer_encoding, titer_scale, titers_base_nolookup):
    # make interquartile range chart

    interquartile_range = (
        titers_base_nolookup.transform_joinaggregate(
            median_titer="median(titer)",
            titer_q1="q1(titer)",
            titer_q3="q3(titer)",
            groupby=["virus"],
        )
        .encode(
            **{titer_encoding: alt.X("titer", scale=titer_scale)},
            tooltip=[
                "virus",
                alt.Tooltip("median_titer:Q", format=".1f"),
                alt.Tooltip("titer_q1:Q", format=".1f"),
                alt.Tooltip("titer_q3:Q", format=".1f"),
                "strain_type:N",
                "subclade:N",
            ],
        )
        .mark_errorband(extent="iqr", opacity=0.5, interpolate="linear")
    )

    # facet_and_add_lookups(interquartile_range + median_points)
    return (interquartile_range,)


@app.cell
def _(
    alt,
    human_sera_plots_params,
    mo,
    titer_encoding,
    titer_lower_limit,
    titers_base_nolookup,
    virus_selection,
):
    # make fraction below titer cutoff chart

    titer_cutoff = human_sera_plots_params["titer_cutoff"]
    mo.output.append(mo.md(f"Setting initial {titer_cutoff=}"))

    titer_cutoff_slider = alt.param(
        value=titer_cutoff,
        bind=alt.binding_range(
            min=titer_lower_limit,
            max=1000,
            step=5,
            name="fraction sera below this cutoff",
        ),
    )

    # make titer cutoff chart
    frac_below_cutoff = (
        titers_base_nolookup.add_params(titer_cutoff_slider)
        .transform_calculate(below_cutoff=alt.datum["titer"] < titer_cutoff_slider)
        .transform_aggregate(
            n_below_cutoff="sum(below_cutoff)",
            n_total="distinct(serum)",
            groupby=["virus", "subtype", "strain_type", "subclade", "group"],
        )
        .transform_calculate(
            frac_below_cutoff=alt.datum["n_below_cutoff"] / alt.datum["n_total"]
        )
        .encode(
            **{
                titer_encoding: alt.X(
                    "frac_below_cutoff:Q", title="fraction below cutoff"
                )
            },
            tooltip=[
                "virus",
                alt.Tooltip("frac_below_cutoff:Q", format=".2f"),
                "strain_type:N",
                "subclade:N",
            ],
            color=alt.condition(virus_selection, alt.value("red"), alt.value("black")),
        )
        .mark_bar(opacity=0.8)
    )

    # facet_and_add_lookups(frac_below_cutoff)
    return frac_below_cutoff, titer_cutoff


@app.cell
def _(
    alt,
    chart_htmls,
    circulating_strain_type,
    color_prop_selection,
    dummy_group_chart,
    facet_and_add_lookups,
    facet_orientation,
    frac_below_cutoff,
    interquartile_range,
    itertools,
    median_points,
    mo,
    serum_lines,
    strain_color_prop,
    viruses,
):
    made_chart = {c: False for c in chart_htmls}

    for _subtype, strain_type, (chart_obj, chart_desc, title) in itertools.product(
        viruses["subtype"].unique(),
        ["recent", "vaccine"],
        [
            (
                serum_lines + median_points,
                "individual_sera",
                "median (points) and per-serum (lines) titers",
            ),
            (
                interquartile_range + median_points,
                "interquartile_range",
                "median (points) and interquartile range titers",
            ),
            (
                frac_below_cutoff,
                "frac_below_cutoff",
                "fraction sera below titer cutoff",
            ),
        ],
    ):
        filesuffix = f"{_subtype}_{strain_type}_{chart_desc}.html"
        filename = [c for c in chart_htmls if filesuffix in c]
        assert (
            len(filename) == 1
        ), f"did not find one filesuffix={filesuffix!r} in chart_htmls={chart_htmls!r}"
        filename = filename[0]

        # strain types to plot
        strain_types = {
            "recent": [circulating_strain_type, "recent_vaccine"],
            "vaccine": ["vaccine", "recent_vaccine"],
        }[strain_type]

        # ---- Make the legend for the colored strain labels ------------------------
        # viruses plotted
        plotted_viruses = (
            viruses[
                (viruses["subtype"] == _subtype)
                & (viruses["strain_type"].isin(strain_types))
            ]["virus"]
            .unique()
            .tolist()
        )
        # get the virus colors plotted for the labels
        plotted_colors = strain_color_prop.query("virus in @plotted_viruses")[
            ["color_prop", "color"]
        ].drop_duplicates()

        label_color_legend = (
            alt.Chart(plotted_colors)
            .add_params(color_prop_selection)
            .mark_point(opacity=0)  # invisible mark; we just want the legend
            .encode(
                fill=alt.Fill(
                    "color_prop",
                    title="virus type (click to select)",
                    scale=alt.Scale(
                        domain=list(reversed(plotted_colors["color_prop"].tolist())),
                        range=list(reversed(plotted_colors["color"].tolist())),
                    ),
                    legend=alt.Legend(symbolType="square"),
                )
            )
            .properties(width=1, height=1)  # tiny plot; legend renders outside
        )
        # ---- Finished making the legend for the colored strain labels ---------------------

        chart = (
            alt.vconcat(
                (
                    facet_and_add_lookups(chart_obj)
                    .transform_filter(alt.datum["subtype"] == _subtype)
                    .transform_filter(
                        alt.FieldOneOfPredicate("strain_type", strain_types)
                    )
                ),
                label_color_legend,
                dummy_group_chart,
                spacing=1,
            )
            .resolve_scale(fill="independent")
            .configure_axis(
                grid=False,
                titleFontWeight="normal",
                titleFontSize=13,
                labelOverlap=True,
            )
            .configure_header(
                title=None,
                labelOrient="top" if facet_orientation == "vertical" else "right",
                labelFontSize=13,
                labelPadding=2,
            )
            .configure_view(stroke="black")
            .configure_facet(spacing=8)
            .configure_legend(
                labelFontSize=12,
                titleFontSize=13,
                symbolStrokeWidth=1,
                symbolOpacity=1,
                symbolStrokeColor="black",
                columns=12,
                orient="bottom",
            )
            .properties(
                title=alt.TitleParams(
                    f"{title} for {_subtype} {strain_type} strains",
                    anchor="middle",
                    fontSize=13,
                )
            )
        )
        if not any(made_chart.values()):
            mo.output.append(
                mo.md("Displaying just the first chart here (since they are large).")
            )
            mo.output.append(chart)

        mo.output.append(mo.md(f"Saving to filename={filename!r}\n"))
        chart.save(filename)

        made_chart[filename] = True

    assert all(made_chart.values()), f"made_chart={made_chart!r}"
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Make data frame / CSV of summarized titers
    The above plots summarize the per-virus titers with several statistics: medians, interquartile range, and fraction below cutoff.
    Here we make and write a data frame with these summarized values.
    """
    )
    return


@app.cell
def _(mo, pd, summarized_titers_csv, titer_cutoff, titers, viruses):
    summarized_titers = (
        pd.concat([titers.assign(group="All"), titers])
        .merge(viruses, on=["virus"], how="left", validate="many_to_one")
        .groupby(
            ["subtype", "strain_type", "subclade", "virus", "group"], as_index=False
        )
        .aggregate(
            median_titer=pd.NamedAgg("titer", "median"),
            titer_q1=pd.NamedAgg("titer", lambda s: s.quantile(0.25)),
            titer_q3=pd.NamedAgg("titer", lambda s: s.quantile(0.75)),
            frac_below_cutoff=pd.NamedAgg(
                "titer", lambda s: (s < titer_cutoff).sum() / len(s)
            ),
        )
        .rename(
            columns={
                "frac_below_cutoff": f"frac_w_titer_below_{titer_cutoff}",
                "group": "serum_group",
            }
        )
    )

    mo.output.append(mo.md(f"Saving summarized titers to {summarized_titers_csv=}"))
    summarized_titers.to_csv(summarized_titers_csv, float_format="%.3f", index=False)

    mo.output.append(summarized_titers)
    return


if __name__ == "__main__":
    app.run()
