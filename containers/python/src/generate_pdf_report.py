import sys
import json
from functools import reduce

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import A3, A4, landscape, letter, portrait
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import (
    Image,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)
from svglib.svglib import svg2rlg

"""
Data import and manipulation
"""




def generate(input, output):

    with open(input, "r") as fp:
        input_data = json.load(fp)



    metadata = input_data["Metadata"]
    varaints_data = input_data["Locus Results"]

    # df = pd.DataFrame(varaints_data).head(20)
    df = pd.DataFrame(varaints_data)

    df["Coverage"] = df["Coverage"].map(lambda x: round(x, 4))

    # Calculate max length of each column
    max_cell_length = df.map(str).map(len).max()
    max_header_length = df.columns.map(str).map(len)
    max_length = reduce(
        lambda acc, x: acc + [max(x[0], x[1])],
        zip(list(max_cell_length), list(max_header_length)),
        [],
    )

    """
    PDF setup

    """

    # Create PDF document
    # PDF_SIZE = landscape(A3)
    # A4 -> 8.3 x 11.7 in
    PDF_SIZE = portrait(A4)
    pdf = SimpleDocTemplate(output, pagesize=PDF_SIZE)
    available_width = PDF_SIZE[0] - 2 * pdf.leftMargin
    content = []

    PRIMARY_COLOR = "#0c4977"
    SECONDARY_COLOR = "#90c85c"
    TERTIARY_COLOR = "#ffffff"
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        "TitleStyle", parent=styles["Title"], textColor=PRIMARY_COLOR
    )
    subtitle_style = ParagraphStyle(
        "SubtitleStyle", parent=styles["Heading2"], textColor=PRIMARY_COLOR
    )
    body_style = styles["BodyText"]

    desired_width = PDF_SIZE[0]  # Target width in points (1 point = 1/72 inch)

    # # Report logo
    # logo_path = "ADL Logo.PNG"
    # report_logo = Image(logo_path)
    # print(report_logo.drawWidth, report_logo.drawHeight)
    # ratio = report_logo.drawWidth / report_logo.drawHeight
    # report_logo.drawWidth = 100
    # report_logo.drawHeight = 100 / ratio
    # # report_logo = Image(logo_path, width=100, height=50)
    # content.append(report_logo)
    # content.append(Spacer(1, 36))

    # Add a title
    content.append(Paragraph("First Report", title_style))
    content.append(Spacer(1, 12))

    # Add a section with subtitle and a paragraph
    content.append(Paragraph("Overview", subtitle_style))
    intro_text = (
        "This document itemizes STR genotypes generated using DRAGEN's ExpansionHunter genotyping tool for the following STR Loci: "
        + ", ".join(
            df["Locus ID"].unique()
        )  # TODO: Check if this is correct sort order
    )
    content.append(Paragraph(intro_text, body_style))
    content.append(Spacer(1, 12))
    content.append(Paragraph("Sample information:", subtitle_style))
    content.append(Paragraph(f"Sample ID: {metadata['Sample ID']}", body_style))
    content.append(Paragraph(f"Sex: {metadata['Sex']}", body_style))
    content.append(Spacer(1, 12))
    content.append(Paragraph("Locus Results:", subtitle_style))
    content.append(Spacer(1, 12))

    """
    Table creation and styling

    """
    WORD_SIZE = 6  # Adjust based on font and font size
    column_widths = [length * WORD_SIZE for length in max_length]
    print(PDF_SIZE)
    total_width = sum(column_widths)
    print(total_width, available_width)
    # Check if total width exceeds available width, adjust if necessary
    if total_width > available_width:
        # Example adjustment: scale down proportionally (simple method)
        scale_factor = available_width / total_width
        column_widths = [width * scale_factor for width in column_widths]

    table_style = TableStyle(
        [
            ("FONT", (0, 0), (-1, -1), "Helvetica", 4),  # Default font for the table
            ("FONT", (1, 0), (1, -1), "Helvetica-Oblique", 4),  # Italic
            ("BACKGROUND", (0, 0), (-1, 0), PRIMARY_COLOR),
            ("TEXTCOLOR", (0, 0), (-1, 0), TERTIARY_COLOR),
            ("ALIGN", (0, 0), (-1, -1), "CENTER"),
            ("FONT", (0, 0), (-1, 0), "Helvetica-Bold", 5),
            # ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("BOTTOMPADDING", (0, 0), (-1, 0), 3),
            ("TOPPADDING", (0, 0), (-1, 0), 3),
            ("BACKGROUND", (0, 1), (-1, -1), TERTIARY_COLOR),
            ("GRID", (0, 0), (-1, -1), 1, TERTIARY_COLOR),
        ]
    )

    # Create table
    tables_scheme = [
        (
            df.sort_values("Locus ID"),
            "Sorted by Locus ID",
        ),
        (
            df.sort_values("Coverage", ascending=False),
            "Sorted by Coverage",
        ),
        (
            df.assign(
                chr_num=df["Reference Region"].str.extract(r"(\d+):").astype("Int64"),
                start_pos=df["Reference Region"].str.extract(r"-(\d+)").astype("Int64"),
            )
            .sort_values(by=["chr_num", "start_pos"])
            .drop(["chr_num", "start_pos"], axis=1),
            "Sorted by Reference Region",
        ),
        (
            df.sort_values("Genotype", ascending=False),
            "Sorted by Genotype",
        ),
    ]

    for i, (table_df, title) in enumerate(tables_scheme):
        # Convert DataFrame to a list of lists for ReportLab
        data_for_report = [table_df.columns.to_list()] + table_df.values.tolist()
        table = Table(data_for_report, colWidths=column_widths, style=table_style)
        # Add table to PDF document
        if i in [1, 3]:
            content.append(PageBreak())
        content.append(Paragraph(title, body_style))
        content.append(Spacer(1, 6))
        content.append(table)
        content.append(Spacer(1, 24))

        # content.append(PageBreak())

    from tempfile import TemporaryDirectory

    from plots import get_fig1, get_fig2, get_fig3, get_fig4, get_fig5, get_fig6

    with TemporaryDirectory() as tmpdir:

        for func in [get_fig1 , get_fig3, get_fig2, get_fig4, get_fig5, get_fig6]:
            plotly_fig = func(df)
            filename = f"{tmpdir}/{func.__name__}.svg"
            plotly_fig.write_image(filename, format="svg")
            fig = svg2rlg(filename)

            scale_factor = available_width / fig.width

            # Apply scaling
            fig.width *= scale_factor
            fig.height *= scale_factor
            fig.scale(scale_factor, scale_factor)
            # print(fig.width, fig.height)
            content.append(fig)

    # Footnote
    footnote_style = ParagraphStyle(
        "FootnoteStyle", parent=styles["Normal"], fontSize=8, textColor=colors.grey
    )
    footnote = Paragraph(
        "CLIA #11D2258259, Lab Director Dr.Hubert Pare.", footnote_style
    )
    content.append(footnote)

    # Build PDF
    pdf.build(content)

    print(f"PDF report '{output}' has been created.")


# if __name__ == "__main__":

#     import sys

#     arguments = sys.argv[1:]
#     # -i, sample_id_raw.json, ...
#     reports = [
#         # (
#         #     "result_from_json_big.json",
#         #     "persida_first_report_big.pdf",
#         # ),
#         (
#             "result_from_json.json",
#             "persida_first_report.pdf",
#         ),
#     ]
#     for report in reports:
#         generate(*report)

def main():
    input = sys.argv[1]
    output = sys.argv[2]
    # print(1/0) PUKLO
    generate(input, output)
    # print(1/0) NIJE PUKLO
    
if __name__ == "__main__":
    main()


# if __name__ == "__main__":

#     import sys

#     # arguments = sys.argv[1:]
#     # -i, sample_id_raw.json, ...
#     reports = [
#         # (
#         #     "result_from_json_big.json",
#         #     "persida_first_report_big.pdf",
#         # ),
#         (
#             sys.argv[1],
#             sys.argv[2],
#         ),
#     ]

#     for report in reports:
#         generate(report)


