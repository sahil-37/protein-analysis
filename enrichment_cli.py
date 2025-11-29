#!/usr/bin/env python3
"""
Unified Protein Enrichment CLI

Generate comprehensive protein enrichment reports with:
- UniProt metadata
- Biophysical properties (MW, pI, GRAVY, etc.)
- Charge profile with pH-dependent visualization
- Structural feature predictions
- JSON and PDF output

Usage:
    python enrichment_cli.py P04637
    python enrichment_cli.py P04637 --pdf
    python enrichment_cli.py P04637 P68871 P15018 --pdf
    python enrichment_cli.py P04637 --json --pdf
"""

import argparse
import sys
from pathlib import Path
from protein_enrichment_service import ProteinEnrichmentService
import json


def generate_interactive_charge_profile(profile, output_filename=None):
    """
    Generate an interactive HTML charge profile with tooltips using Plotly

    Args:
        profile: ProteinEnrichmentProfile object
        output_filename: Output HTML filename

    Returns:
        Path to generated HTML file
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("  ⚠️  Plotly not installed - skipping interactive HTML. Run: pip install plotly")
        return None

    if output_filename is None:
        output_filename = f"{profile.general.uniprot_id}_charge_profile_interactive.html"

    ph_vals = profile.charges.ph_values
    charges = profile.charges.charges
    pi = profile.charges.isoelectric_point
    charge_74 = profile.charges.charge_at_ph_7_4

    # Create interactive plot
    fig = go.Figure()

    # Add charge profile curve with custom hover text
    hover_text = [
        f"<b>pH: {ph:.1f}</b><br>Net Charge: {charge:+.2f} e⁻"
        for ph, charge in zip(ph_vals, charges)
    ]

    fig.add_trace(go.Scatter(
        x=ph_vals,
        y=charges,
        mode='lines',
        name='Net Charge',
        line=dict(color='#2E86AB', width=3),
        hovertext=hover_text,
        hoverinfo='text',
        fill='tozeroy',
        fillcolor='rgba(46, 134, 171, 0.15)',
    ))

    # Add pI line
    fig.add_vline(
        x=pi,
        line_dash="dash",
        line_color="#E63946",
        line_width=2,
        annotation_text=f"pI = {pi:.2f}",
        annotation_position="top",
    )

    # Add physiological pH line
    fig.add_vline(
        x=7.4,
        line_dash="dot",
        line_color="#06A77D",
        line_width=2,
        annotation_text=f"pH 7.4 (Charge: {charge_74:+.2f})",
        annotation_position="top",
    )

    # Add zero charge line
    fig.add_hline(
        y=0,
        line_color="#999999",
        line_width=1,
        line_dash="solid",
        opacity=0.5,
    )

    # Update layout
    fig.update_layout(
        title=f"Interactive Charge Profile - {profile.general.uniprot_id} ({profile.general.protein_name})",
        xaxis_title="pH",
        yaxis_title="Net Charge (e⁻)",
        hovermode='x unified',
        template='plotly_white',
        height=600,
        width=1000,
        font=dict(size=12),
        showlegend=True,
    )

    # Save to HTML
    fig.write_html(output_filename)
    return output_filename


def generate_charge_profile_pdf(profile, output_filename=None):
    """
    Generate a standalone PDF report with charge profile graph

    Args:
        profile: ProteinEnrichmentProfile object
        output_filename: Output PDF filename

    Returns:
        Path to generated PDF
    """
    from protein_design_pdf_export import create_charge_profile_chart
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.units import inch
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
    from datetime import datetime

    if output_filename is None:
        output_filename = f"{profile.general.uniprot_id}_enrichment_report.pdf"

    # Create PDF document
    doc = SimpleDocTemplate(
        output_filename,
        pagesize=letter,
        rightMargin=0.75*inch,
        leftMargin=0.75*inch,
        topMargin=0.75*inch,
        bottomMargin=0.75*inch
    )

    # Setup styles
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        name='CustomTitle',
        parent=styles['Heading1'],
        fontSize=24,
        textColor=colors.HexColor('#1F4788'),
        spaceAfter=30,
        alignment=TA_CENTER,
        fontName='Helvetica-Bold'
    )

    section_style = ParagraphStyle(
        name='SectionHeading',
        parent=styles['Heading2'],
        fontSize=14,
        textColor=colors.HexColor('#2E86AB'),
        spaceAfter=12,
        spaceBefore=12,
        fontName='Helvetica-Bold'
    )

    body_style = ParagraphStyle(
        name='CustomBody',
        parent=styles['BodyText'],
        fontSize=10,
        alignment=TA_JUSTIFY,
        spaceAfter=10,
        leading=14
    )

    # Build document content
    elements = []

    # Title
    elements.append(Spacer(1, 0.3*inch))
    elements.append(Paragraph("Protein Enrichment Report", title_style))
    elements.append(Spacer(1, 0.25*inch))

    # Protein Info
    info_data = [
        [f"<b>Protein:</b>", f"{profile.general.protein_name}"],
        [f"<b>Gene:</b>", f"{profile.general.gene_name or 'N/A'}"],
        [f"<b>Accession ID:</b>", f"{profile.general.uniprot_id}"],
        [f"<b>Organism:</b>", f"{profile.general.organism}"],
        [f"<b>Sequence Length:</b>", f"{profile.general.sequence_length} aa"],
        [f"<b>Status:</b>", "Swiss-Prot (Reviewed)" if profile.general.reviewed else "TrEMBL (Unreviewed)"],
        [f"<b>Generated:</b>", datetime.now().strftime("%Y-%m-%d %H:%M:%S")],
    ]

    info_table = Table(info_data, colWidths=[1.8*inch, 3.8*inch])
    info_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (0, -1), colors.HexColor('#E8F4F8')),
        ('ALIGN', (0, 0), (0, -1), 'RIGHT'),
        ('ALIGN', (1, 0), (1, -1), 'LEFT'),
        ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, -1), 10),
        ('GRID', (0, 0), (-1, -1), 1, colors.HexColor('#BDC3C7')),
        ('ROWBACKGROUNDS', (0, 0), (-1, -1), [colors.white, colors.HexColor('#F0F5F7')]),
        ('TOPPADDING', (0, 0), (-1, -1), 8),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
    ]))

    elements.append(info_table)
    elements.append(Spacer(1, 0.3*inch))

    # Biophysical Properties
    elements.append(Paragraph("Biophysical Properties", section_style))

    biophys = profile.biophysical
    biophys_data = [
        ['Property', 'Value', 'Interpretation'],
        ['Molecular Weight', f"{biophys.molecular_weight_kda:.2f} kDa", 'Protein mass'],
        ['Isoelectric Point (pI)', f"{biophys.isoelectric_point_pi:.2f}", 'pH of zero charge'],
        ['Aromaticity', f"{biophys.aromaticity*100:.2f}%", f"{(biophys.aromaticity)*100:.1f}% aromatic residues (F,W,Y)"],
        ['GRAVY Index', f"{biophys.gravy:+.3f}", "Hydropathy (negative=hydrophilic)"],
        ['Instability Index', f"{biophys.instability_index:.1f}", "Stability prediction (>40=unstable)"],
        ['Aliphatic Index', f"{biophys.aliphatic_index:.1f}", "Thermostability index"],
        ['Cysteine Count', f"{biophys.cysteine_count}", f"Total cysteines in sequence"],
    ]

    biophys_table = Table(biophys_data, colWidths=[2.0*inch, 1.2*inch, 2.3*inch])
    biophys_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#2E86AB')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 11),
        ('FONTSIZE', (0, 1), (-1, -1), 9),
        ('BACKGROUND', (0, 1), (-1, -1), colors.HexColor('#F0F5F7')),
        ('GRID', (0, 0), (-1, -1), 1, colors.HexColor('#BDC3C7')),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F0F5F7')]),
        ('TOPPADDING', (0, 0), (-1, -1), 8),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
    ]))

    elements.append(biophys_table)
    elements.append(Spacer(1, 0.3*inch))

    # Disulfide Bonds (from UniProt)
    if profile.predicted_features and profile.predicted_features.disulfide_bonds:
        elements.append(Paragraph("Confirmed Disulfide Bonds", section_style))

        disulfide_data = [
            ['Bond #', 'Cys Position 1', 'Cys Position 2', 'Source'],
        ]

        for i, bond in enumerate(profile.predicted_features.disulfide_bonds, 1):
            source = "UniProt (Experimentally Confirmed)" if bond.is_experimentally_confirmed else "Predicted"
            disulfide_data.append([
                str(i),
                f"Cys{bond.cys_position_1}",
                f"Cys{bond.cys_position_2}",
                source
            ])

        disulfide_table = Table(disulfide_data, colWidths=[0.8*inch, 1.2*inch, 1.2*inch, 2.3*inch])
        disulfide_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#2E86AB')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('ALIGN', (3, 0), (3, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 11),
            ('FONTSIZE', (0, 1), (-1, -1), 9),
            ('BACKGROUND', (0, 1), (-1, -1), colors.HexColor('#F0F5F7')),
            ('GRID', (0, 0), (-1, -1), 1, colors.HexColor('#BDC3C7')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F0F5F7')]),
            ('TOPPADDING', (0, 0), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
        ]))

        elements.append(disulfide_table)
        elements.append(Spacer(1, 0.3*inch))

    # Charge Profile
    elements.append(Paragraph("pH-Dependent Charge Profile", section_style))
    elements.append(Spacer(1, 0.15*inch))

    # Create and add charge profile chart
    ph_vals = profile.charges.ph_values
    charges = profile.charges.charges
    pi = profile.charges.isoelectric_point
    charge_74 = profile.charges.charge_at_ph_7_4

    chart = create_charge_profile_chart(ph_vals, charges, pi)
    elements.append(chart)
    elements.append(Spacer(1, 0.2*inch))

    # Charge interpretation (facts only - no generic applications)
    charge_text = f"""
    The protein has an isoelectric point (pI) of <b>{pi:.2f}</b>, meaning it carries zero net charge at this pH.
    At physiological pH (7.4), the protein carries a net charge of <b>{charge_74:+.2f} e⁻</b>.<br/>
    <br/>
    <b>Charge characteristics:</b><br/>
    • Below pI (pH &lt; {pi:.2f}): Protein is positively charged<br/>
    • Above pI (pH &gt; {pi:.2f}): Protein is negatively charged<br/>
    • At pH 7.4: Charge is {charge_74:+.2f} ({"negative" if charge_74 < 0 else "positive" if charge_74 > 0 else "neutral"})<br/>
    """
    elements.append(Paragraph(charge_text, body_style))

    # Build PDF
    doc.build(elements)

    return output_filename


def main():
    """Main CLI entry point"""

    parser = argparse.ArgumentParser(
        description='Unified Protein Enrichment Service - Generate comprehensive protein analysis reports',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Terminal output only
  python enrichment_cli.py P04637

  # Generate JSON and PDF
  python enrichment_cli.py P04637 --json --pdf

  # Batch process multiple proteins
  python enrichment_cli.py P04637 P68871 P15018 --pdf

  # Custom output filename
  python enrichment_cli.py P04637 --pdf -o my_p53_report.pdf
        """
    )

    parser.add_argument('accession_ids', nargs='+',
                       help='UniProt accession ID(s) to enrich')

    parser.add_argument('--json', action='store_true',
                       help='Export results to JSON format')

    parser.add_argument('--pdf', action='store_true',
                       help='Generate PDF report with charge profile graph')

    parser.add_argument('-o', '--output', type=str,
                       help='Output filename (for single protein)')

    args = parser.parse_args()

    # Initialize service
    service = ProteinEnrichmentService()

    # Process each protein
    for acc_id in args.accession_ids:
        print(f"\n{'='*60}")
        print(f"Processing: {acc_id}")
        print('='*60)

        try:
            # Enrich protein
            profile = service.enrich_protein(acc_id)

            # Terminal output
            print(f"\n📊 ENRICHMENT PROFILE: {acc_id}")
            print(f"   Protein: {profile.general.protein_name}")
            print(f"   Length: {profile.general.sequence_length} aa")
            print(f"   MW: {profile.biophysical.molecular_weight_kda:.2f} kDa")
            print(f"   pI: {profile.biophysical.isoelectric_point_pi:.2f}")
            print(f"   Charge @ pH 7.4: {profile.charges.charge_at_ph_7_4:+.2f}")
            print(f"   Aromaticity: {profile.biophysical.aromaticity*100:.1f}%")
            print(f"   GRAVY: {profile.biophysical.gravy:+.3f}")
            print(f"   Instability Index: {profile.biophysical.instability_index:.1f}")

            # JSON export
            if args.json:
                json_file = args.output or f"{acc_id}_enriched.json"
                profile.to_json(json_file)
                print(f"   ✓ JSON exported: {json_file}")

            # PDF export
            if args.pdf:
                pdf_file = args.output or f"{acc_id}_enrichment_report.pdf"
                pdf_file = generate_charge_profile_pdf(profile, pdf_file)
                print(f"   ✓ PDF generated: {pdf_file}")

                # Also generate interactive HTML charge profile
                html_file = generate_interactive_charge_profile(profile)
                if html_file:
                    print(f"   ✓ Interactive HTML generated: {html_file}")

        except Exception as e:
            print(f"   ✗ Error: {e}")
            import traceback
            traceback.print_exc()

    print(f"\n{'='*60}")
    print("✓ Enrichment complete")
    print('='*60)


if __name__ == '__main__':
    main()
