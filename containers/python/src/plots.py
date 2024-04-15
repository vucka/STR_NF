import json

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Dash, dcc, html

"""
axis -> ['Locus ID', 'Allele Count', 'Coverage', 'Fragment Length', 'Read Length', 'Reference Region', 'Variant Type', 'Repeat Unit', 'Genotype']
barmode -> ['stack', 'group', 'overlay', 'relative']
"""


PRIMARY_COLOR = "#0c4977"
SECONDARY_COLOR = "#90c85c"
TERTIARY_COLOR = "#ffffff"

COLORS_TO_USE_DISCRETE = px.colors.diverging.balance
COLORS_TO_USE_CONTINUOUS = px.colors.sequential.Blues


def modify_df(func):
    def wrap(df, *args, **kwargs):
        df = df.assign(
            Chromosome=df["Reference Region"].str.extract(r"(\d+):").astype("Int64"),
            Chromosome_str=df["Reference Region"].str.extract(r"(chr\d+):"),
            Position=df["Reference Region"].str.extract(r"-(\d+)").astype("Int64"),
            Start_End_Position=df["Reference Region"].str.extract(r"(\d+-\d+)"),
            Locus_ID_Repeat=df["Locus ID"].str.cat(df["Repeat Unit"], sep=": "),
            Genotype_repeat_allele_1=df["Genotype"]
            .str.extract(r"(\d+)/")
            .astype("Int64"),
            Genotype_repeat_allele_2=df["Genotype"]
            .str.extract(r"/(\d+)")
            .astype("Int64"),
            confidence_interval_min_1=df["Genotype Confidence Interval"]
            .str.extract(r"(\d+)-")
            .astype("Int64"),
            confidence_interval_max_1=df["Genotype Confidence Interval"]
            .str.extract(r"-(\d+)/")
            .astype("Int64"),
            confidence_interval_min_2=df["Genotype Confidence Interval"]
            .str.extract(r"/(\d+)-")
            .astype("Int64"),
            confidence_interval_max_2=df["Genotype Confidence Interval"]
            .str.extract(r"-(\d+)$")
            .astype("Int64"),
        )
        df = df.assign(
            Locus_ID_Chromosome=df["Locus ID"].str.cat(df["Chromosome_str"], sep=": "),
            Locus_ID_Repeat_Genotype=df["Locus_ID_Repeat"].str.cat(
                df["Genotype"], sep=" - "
            ),
        )
        df = df.sort_values(["Coverage"], ascending=False)
        # print(df.head(20))
        # df = df.head(20)
        # df_chrom_sorted = df.sort_values(by=["Chromosome", "Position"])
        return func(df, *args, **kwargs)

    return wrap


@modify_df
def get_fig1(df, colors=COLORS_TO_USE_CONTINUOUS):
    fig = px.bar(
        df,
        x="Locus_ID_Repeat",
        y="Coverage",
        color="Fragment Length",
        text="Chromosome_str",
        color_continuous_scale=colors,
    )
    fig.update_traces(
        marker_line_color="black",
        marker_line_width=0.5,
        opacity=0.6,
        textfont_size=12,
        textangle=0,
        textposition="outside",
        cliponaxis=False,
        showlegend=False,
    )
    fig.update_layout(
        xaxis_title="",
        yaxis_title="",
        title=dict(
            text="Coverage by Locus ID: Repeat Unit (Color: Fragment Length, Text: Chromosome)",
            # font=dict(size=25),
        ),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )
    fig.add_annotation(
        text="Locus ID: Repeat Unit",
        x=0.5,
        y=1.05,
        showarrow=False,
        xref="paper",
        yref="paper",
        font=dict(
            size=20,
        ),
    )
    fig.add_annotation(
        text="Coverage",
        x=-0.1,
        y=0.5,
        showarrow=False,
        xref="paper",
        yref="paper",
        font=dict(
            size=20,
        ),
        textangle=-90,
    )

    fig.update(layout_coloraxis_showscale=False)
    return fig


@modify_df
def get_fig2(df, colors=COLORS_TO_USE_CONTINUOUS):
    PADDING = 15
    df.sort_values(by=["Coverage"], inplace=True)
    minimum_coverage = max(0, df["Coverage"].min() - PADDING)
    if minimum_coverage:
        text_position = "inside"
    else:
        text_position = "outside"

    fig = px.bar(
        df,
        x="Coverage",
        y="Locus_ID_Repeat",
        color="Fragment Length",
        text="Reference Region",
        color_continuous_scale=colors,
    )
    fig.update_traces(
        marker_line_color="black",
        marker_line_width=0.5,
        opacity=0.6,
        textfont_size=12,
        textangle=0,
        textposition=text_position,
        cliponaxis=False,
    )
    fig.update_layout(
        title="Locus ID: Repeat Unit by Coverage (Color: Fragment Length, Text: Reference Region)",
        # height=800,
        xaxis_title="",
        yaxis_title="",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )
    fig.update_coloraxes(
        colorbar=dict(
            title="Fragment Length",
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1,
        )
    )
    fig.add_annotation(
        text="Coverage",
        x=0.5,
        y=-0.1,
        showarrow=False,
        xref="paper",
        yref="paper",
        font=dict(
            size=20,
        ),
    )
    fig.add_annotation(
        text="Locus ID: Repeat Unit",
        x=-0.3,
        y=0.5,
        showarrow=False,
        xref="paper",
        yref="paper",
        font=dict(
            size=20,
        ),
        textangle=-90,
    )
    fig.update(layout_coloraxis_showscale=False)

    return fig


@modify_df
def get_fig3(df, colors=COLORS_TO_USE_DISCRETE):
    fig = px.scatter(
        df,
        x="Fragment Length",
        y="Coverage",
        color="Locus_ID_Repeat_Genotype",
        size=df["Repeat Unit"].apply(lambda x: 0.001 * (len(x) ** 2)),
        # text="Locus_ID_Repeat",
        color_discrete_sequence=colors,
    )
    fig.update_traces(textposition="top center", cliponaxis=False)
    # fig.update_traces(textfont_size=12, textposition="top center", cliponaxis=False)
    fig.update_layout(
        title="Coverage by Fragment Length (Color: Genotype, Size: Repeat Unit)",
        xaxis_title="Fragment Length",
        yaxis_title="Coverage",
        legend_title="Locus ID: Repeat Unit - Genotype",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )
    return fig


@modify_df
def get_fig4(df, color_scale=COLORS_TO_USE_DISCRETE):

    repeat_units = list(df["Repeat Unit"].unique())
    # combine Locus ID and Genotype to one string for labels

    from collections import defaultdict

    r_colors = ["/"] * len(repeat_units)
    l_colors = df["Genotype"].to_list()

    col_dict = {}
    counter = 0
    for i, genotype in enumerate(reversed(r_colors + l_colors)):
        if genotype not in col_dict:
            col_dict[genotype] = counter
            counter += 1

    colors = [col_dict[genotype] for genotype in r_colors + l_colors]
    df["Locus_ID_Chromosome"] = (
        df["Locus_ID_Chromosome"]
        + "<br>Genotype= "
        + df["Genotype"]
        + "<br>Coverage= "
        + round(df["Coverage"], 4).astype(str)
    )

    labels = repeat_units + df["Locus_ID_Chromosome"].to_list()
    parents = [""] * len(repeat_units) + df["Repeat Unit"].to_list()
    c_dict = {k: 0 for k in labels}
    # group by repeat unit
    for i, row in df.iterrows():
        c_dict[row["Locus_ID_Chromosome"]] = row["Coverage"]
    fig = go.Figure(
        data=go.Sunburst(
            labels=labels,
            parents=parents,
            values=list(c_dict.values()),
            marker=dict(
                line=dict(
                    color=PRIMARY_COLOR,
                    width=0.5,
                ),
                colors=colors,
                colorscale="Blues",
            ),
        ),
        layout=go.Layout(
            title="Coverage by Repeat Unit (Color: Genotype)",
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
            width=800,
            height=800,
        ),
    )

    return fig


@modify_df
def get_fig5(df, colors=COLORS_TO_USE_DISCRETE):
    new_df = (
        df.groupby(
            [
                "Locus_ID_Repeat",
                "Genotype_repeat_allele_1",
                "Genotype_repeat_allele_2",
            ]
        )
        .size()
        .reset_index(name="Count")
    )
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=new_df["Locus_ID_Repeat"],
            y=new_df["Genotype_repeat_allele_1"],
            name="first allele",
            text=new_df["Genotype_repeat_allele_1"],
            marker_color=colors[1],
        )
    )
    fig.add_trace(
        go.Bar(
            x=new_df["Locus_ID_Repeat"],
            y=new_df["Genotype_repeat_allele_2"],
            name="second allele",
            text=new_df["Genotype_repeat_allele_2"],
            marker_color=colors[3],
        )
    )

    fig.update_traces(
        texttemplate="%{text}",
        textposition="outside",
        cliponaxis=False,
        # textangle=0,
        textfont_size=12,
    )
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(
            title="Locus ID: Repeat ",
            titlefont_size=16,
            tickfont_size=14,
        ),
        yaxis=dict(
            title="Allele Count",
            titlefont_size=16,
            tickfont_size=14,
        ),
    )
    return fig


@modify_df
def get_fig6(df, colors=COLORS_TO_USE_DISCRETE):
    new_df = (
        df.groupby(
            [
                "Locus_ID_Repeat",
                "confidence_interval_min_1",
                "confidence_interval_max_1",
                "confidence_interval_min_2",
                "confidence_interval_max_2",
                "Genotype_repeat_allele_1",
                "Genotype_repeat_allele_2",
            ]
        )
        .size()
        .reset_index(name="Count")
    )
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=new_df["Locus_ID_Repeat"],
            y=new_df["confidence_interval_max_1"] - new_df["confidence_interval_min_1"],
            name="first allele",
            text="Confidence Interval:<br>"
            + new_df["confidence_interval_min_1"].astype(str)
            + " - "
            + new_df["confidence_interval_max_1"].astype(str)
            + "<br> Allele Count: "
            + new_df["Genotype_repeat_allele_1"].astype(str),
            marker_color=colors[1],
            base=new_df["confidence_interval_min_1"],
        )
    )
    fig.add_trace(
        go.Bar(
            x=new_df["Locus_ID_Repeat"],
            y=new_df["confidence_interval_max_2"] - new_df["confidence_interval_min_2"],
            name="second allele",
            text="Confidence Interval:<br>"
            + new_df["confidence_interval_min_2"].astype(str)
            + " - "
            + new_df["confidence_interval_max_2"].astype(str)
            + "<br> Allele Count: "
            + new_df["Genotype_repeat_allele_2"].astype(str),
            marker_color=colors[3],
            base=new_df["confidence_interval_min_2"],
        )
    )
    # fig.add_trace(
    #     go.Bar(
    #         x=new_df["Locus_ID_Repeat"],
    #         y=new_df["Genotype_repeat_allele_2"],
    #         name="second allele",
    #         text=new_df["Genotype_repeat_allele_2"],
    #         marker_color=colors[3],
    #     )
    # )
    # fig.add_trace(
    #     go.Scatter(
    #         x=new_df["Locus_ID_Repeat"],
    #         y=new_df["Genotype_repeat_allele_1"],
    #         mode="markers",
    #         name="first allele",
    #     )
    # )
    # fig.add_trace(
    #     go.Scatter(
    #         x=new_df["Locus_ID_Repeat"],
    #         y=new_df["Genotype_repeat_allele_2"],
    #         mode="markers",
    #         name="second allele",
    #     )
    # )

    fig.update_traces(
        texttemplate="%{text}",
        textposition="outside",
        cliponaxis=False,
        textangle=0,
        textfont_size=12,
    )
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(
            title="Locus ID: Repeat ",
            titlefont_size=16,
            tickfont_size=14,
        ),
        yaxis=dict(
            title="Allele Count Confidence Interval",
            titlefont_size=16,
            tickfont_size=14,
        ),
    )
    return fig


@modify_df
def get_fig7(df, colors=COLORS_TO_USE_DISCRETE):

    new_df = (
        df.groupby(
            [
                "Locus_ID_Repeat",
                "confidence_interval_min_1",
                "confidence_interval_max_1",
                "confidence_interval_min_2",
                "confidence_interval_max_2",
            ]
        )
        .size()
        .reset_index(name="Count")
    )
    fig = go.Figure()
    width = 0.5
    gap = 0.5
    internal_gap = width + 0.1
    interval = 0
    for i, row in new_df.iterrows():
        fig.add_trace(
            go.Scatter(
                x=[
                    i * (width + gap) + interval * internal_gap,
                    i * (width + gap) + interval * internal_gap,
                    i * (width + gap) + width + +interval * internal_gap,
                    i * (width + gap) + width + +interval * internal_gap,
                    i * (width + gap) + +interval * internal_gap,
                ],
                y=[
                    row["confidence_interval_min_1"],
                    row["confidence_interval_max_1"],
                    row["confidence_interval_max_1"],
                    row["confidence_interval_min_1"],
                    row["confidence_interval_min_1"],
                ],
                fill="toself",
            )
        )
        interval += 1
        fig.add_trace(
            go.Scatter(
                x=[
                    i * (width + gap) + interval * internal_gap,
                    i * (width + gap) + interval * internal_gap,
                    i * (width + gap) + width + interval * internal_gap,
                    i * (width + gap) + width + interval * internal_gap,
                    i * (width + gap) + interval * internal_gap,
                ],
                y=[
                    row["confidence_interval_min_2"],
                    row["confidence_interval_max_2"],
                    row["confidence_interval_max_2"],
                    row["confidence_interval_min_2"],
                    row["confidence_interval_min_2"],
                ],
                fill="toself",
            )
        )
        interval += 1
        # fig.add_shape(
        #     type="rect",
        #     x0=i*(width+gap),
        #     x1=i*(width+gap) + width,
        #     y0=row["confidence_interval_min_1"],
        #     y1=row["confidence_interval_max_1"],
        #     # line=dict(color="RoyalBlue"),
        #     label=dict(text=row["Locus_ID_Repeat"]),
        # )

    return fig
    print(new_df)
    # fig = px.bar(
    #     df,
    #     x="Locus_ID_Repeat",
    #     y="Coverage",
    #     color="Genotype",
    #     text="Repeat Unit",
    #     color_discrete_sequence=colors,

    # )


def main():
    app = Dash(__name__)

    with open("result_from_json_big.json", "r") as fp:
        # with open("result_from_json.json", "r") as fp:
        input_data = json.load(fp)

    varaints_data = input_data["Locus Results"]
    df = pd.DataFrame(varaints_data)
    df = df.head(20)
    # fig1 = get_fig1(df)
    # fig2 = get_fig2(df)
    # fig3 = get_fig3(df)
    # fig4 = get_fig4(df)
    fig5 = get_fig5(df)
    fig6 = get_fig6(df)
    app.layout = html.Div(
        [
            # dcc.Graph(figure=fig1),
            # html.Hr(),
            # dcc.Graph(figure=fig2),
            # html.Hr(),
            # dcc.Graph(figure=fig3),
            # html.Hr(),
            # dcc.Graph(figure=fig4),
            # html.Hr(),
            dcc.Graph(figure=fig5),
            html.Hr(),
            dcc.Graph(figure=fig6),
            # html.Hr(),
            # dcc.Graph(figure=fig7),
        ]
    )
    app.run_server(debug=True, port=8050, host="localhost")


if __name__ == "__main__":
    main()
