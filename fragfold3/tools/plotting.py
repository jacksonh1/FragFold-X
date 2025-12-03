import matplotlib.pyplot as plt
import pandas as pd

plt.style.use("fragfold3.local_data.fragfold3_plotstyle")
import matplotlib.pyplot as plt
import matplotlib.axes
import seaborn as sns
import plotly.subplots as sp
import plotly.graph_objs as go
from Bio import PDB
import tempfile
import py3Dmol
import os
import fragfold3.tools.pdb_tools as pdb_tools


def build_mosaic_z_score_plot(figsize=(15, 5)):
    mos_vector = []
    mos_vector.append(["background"] + ["scores"] * 2)
    mos_vector.append(["background"] + ["logo"] * 2)
    fig, axd = plt.subplot_mosaic(mos_vector, figsize=figsize, layout="constrained")
    # plt.tight_layout()
    return fig, axd


def plot_score_bar_plot(ax: matplotlib.axes.Axes, score_list: list, query_seq: str):
    ax.bar(
        list(range(len(score_list))),
        score_list,
    )
    ax = _format_bar_plot(ax, query_seq)
    return ax


def _format_bar_plot(ax, xlabel_sequence: str):
    """format bar plot"""
    _ = ax.set_xticks(
        list(range(len(xlabel_sequence))),
        labels=list(xlabel_sequence),
    )
    ax.set_xlim(-0.5, len(xlabel_sequence) - 0.5)
    return ax


def plotly_fragfold_results(
    csv_file,
    fragment_source_label,
    receptor_labels,
    cols=None,
    xcol="fragment_center",
    weight_set: str | None = None,
):
    if cols is None:
        cols = ["n_contacts", "iptm", "weighted_contacts"]
    df = pd.read_csv(csv_file)
    # get the mean and std of the scores for each fragment center and weights
    temp = (
        df.groupby([xcol, "weights"])
        .agg(
            {
                "n_contacts": ["mean", "std"],
                "iptm": ["mean", "std"],
                "weighted_contacts": ["mean", "std"],
            }
        )
        .reset_index()
    )
    # columns are multi-indexed, so we need to flatten them
    temp.columns = ["-".join(col).strip() for col in temp.columns.values]
    temp = temp.rename(columns={f"{xcol}-": xcol, "weights-": "weights"})
    fig = sp.make_subplots(
        rows=len(cols),
        cols=1,
        subplot_titles=[
            f"fragment source: {fragment_source_label} - receptor(s): {' + '.join(receptor_labels)}"
            for col in cols
        ],
        vertical_spacing=0.15,
    )
    if weight_set is not None:
        weights = [weight_set]
    else:
        weights = temp["weights"].unique()
    colors = [
        ["rgba(218,124,48,100)", "rgba(218,124,48,0.3)"],
        ["rgba(57,106,177,100)", "rgba(57,106,177,0.3)"],
    ]
    for i, col in enumerate(cols):
        color_count = 0
        for weight in weights:
            temp2 = temp[temp["weights"] == weight]
            x = temp2[xcol]
            y = temp2[f"{col}-mean"]
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    name=f"{weight} - {col} mean",
                    line=dict(color=colors[color_count][0], width=2),
                ),
                row=i + 1,
                col=1,
            )
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y + temp2[f"{col}-std"],
                    mode="lines",
                    marker=dict(color="#444"),
                    line=dict(width=0),
                    showlegend=False,
                    name="+1 std",
                ),
                row=i + 1,
                col=1,
            )
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y - temp2[f"{col}-std"],
                    mode="lines",
                    marker=dict(color="#444"),
                    line=dict(width=0),
                    fillcolor=colors[color_count][1],
                    fill="tonexty",
                    showlegend=False,
                    name="-1 std",
                ),
                row=i + 1,
                col=1,
            )
            color_count += 1
        fig.update_yaxes(title_text=col, row=i + 1, col=1)
        fig.update_xaxes(title_text=xcol, row=i + 1, col=1)
    # change background to white with gray gridlines
    fig.update_xaxes(showgrid=True, gridcolor="lightgray")
    fig.update_yaxes(showgrid=True, gridcolor="lightgray")
    fig.update_layout(
        height=800,
        width=800,
        plot_bgcolor="rgba(255, 255, 255, 1)",
        paper_bgcolor="rgba(255, 255, 255, 1)",
        font=dict(size=14),
    )
    # Make layout more compact
    fig.update_layout(margin=dict(l=50, r=50, t=30, b=30), hovermode="x")
    # add an outline around each plot
    for i in range(len(cols)):
        fig["layout"][f"yaxis{i+1}"].update(
            showline=True, linewidth=2, linecolor="black"
        )
        fig["layout"][f"xaxis{i+1}"].update(
            showline=True, linewidth=2, linecolor="black"
        )
        # add the x-axis ticks
        # fig['layout'][f'xaxis{i+1}'].update(tickmode='linear', dtick=5)
    # reorder the trace labels so that they are in the same order as the plots
    fig.for_each_trace(lambda t: t.update(name=t.name.split(" - ")[0]))
    return fig, df


def plot_fragfold_results(df, title=None, cols=None, xcol="fragment_center", weight_set=None):
    if cols is None:
        cols = ["n_contacts", "iptm", "weighted_contacts"]
    fig, axes = plt.subplots(nrows=3, figsize=(10, 10))
    if weight_set is not None:
        df = df[df["weights"] == weight_set]
    for col, ax in zip(cols, axes.flatten()):
        if "weights" in df.columns:
            sns.lineplot(data=df, x=xcol, y=col, hue="weights", ax=ax)
        else:
            sns.lineplot(data=df, x=xcol, y=col, ax=ax)
        ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
        if title is not None:
            ax.set_title(title)
    fig.tight_layout()
    return fig, axes


def highlight_regions(axes, start, end, **kwargs):
    for ax in axes:
        ax.axvspan(start, end, color="red", alpha=0.3, **kwargs)
    return axes


def biopython_Structure_to_string(structure):
    """Convert a Biopython Structure object to a PDB string."""
    io = PDB.PDBIO()
    io.set_structure(structure)
    temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
    io.save(temp_file.name)
    temp_file.close()
    with open(temp_file.name, "r") as f:
        pdb_string = f.read()
    os.remove(temp_file.name)  # Clean up the temporary file
    return pdb_string


def visualize_pdb_structures(pdb_files):
    # Load structures
    structures = []
    for file in pdb_files:
        parser = PDB.PDBParser()
        structure = parser.get_structure("Structure", file)
        structures.append(structure)

    # Align structures if necessary (example)
    ref_structure = structures[0]
    super_imposer = PDB.Superimposer()

    for structure in structures[1:]:
        atoms_ref = list(ref_structure.get_atoms())
        atoms_structure = list(structure.get_atoms())

        # Ensure both lists have the same length for alignment
        min_len = min(len(atoms_ref), len(atoms_structure))
        atoms_ref = atoms_ref[:min_len]
        atoms_structure = atoms_structure[:min_len]

        super_imposer.set_atoms(atoms_ref, atoms_structure)
        super_imposer.apply(structure)

    view = py3Dmol.view(height=600, width=800)
    for structure in structures:
        pdb_string = biopython_Structure_to_string(structure)
        view.addModel(pdb_string, "pdb")
    view.setStyle({"cartoon": {"colorscheme": "chain"}})
    view.zoomTo()
    return view


def crude_pdb_structure_alignment(pdb_files, chains_to_align):
    structures = []
    for file in pdb_files:
        parser = PDB.PDBParser()
        structure = parser.get_structure("Structure", file)
        structures.append(structure)

    # Align structures if necessary (example)
    ref_structure = structures[0]
    ref_chain_id = chains_to_align[0]
    super_imposer = PDB.Superimposer()
    for structure, chain_id in zip(structures[1:], chains_to_align[1:]):
        atoms_ref = [
            atom
            for atom in ref_structure.get_atoms()
            if atom.get_parent().get_parent().get_id() == ref_chain_id
        ]
        atoms_structure = [
            atom
            for atom in structure.get_atoms()
            if atom.get_parent().get_parent().get_id() == chain_id
        ]
        # Ensure both lists have the same length for alignment
        min_len = min(len(atoms_ref), len(atoms_structure))
        atoms_ref = atoms_ref[:min_len]
        atoms_structure = atoms_structure[:min_len]
        super_imposer.set_atoms(atoms_ref, atoms_structure)
        super_imposer.apply(structure)
    return structures


def visualize_pdb_structures_chain_aln(
    pdb_files, chains_to_align=None, align=True, opacity=1.0
):
    if chains_to_align is None:
        # chains_to_align = ["A"] * len(
        #     pdb_files
        # )  # Default to 'A' for all if not provided
        chains_to_align = [
            pdb_tools.get_chains_from_structure(pdb_files[0])[:-1]
        ] * len(pdb_files)
        # print(chains_to_align)
    if isinstance(chains_to_align, str):
        chains_to_align = [chains_to_align] * len(pdb_files)
    if len(chains_to_align) != len(pdb_files):
        raise ValueError(
            "Length of chains_to_align must match the number of pdb_files."
        )
    if not align:
        structures = []
        for file in pdb_files:
            parser = PDB.PDBParser()
            structure = parser.get_structure("Structure", file)
            structures.append(structure)
    else:
        structures = crude_pdb_structure_alignment(pdb_files, chains_to_align)
    view = py3Dmol.view(height=600, width=800)
    for structure in structures:
        pdb_string = biopython_Structure_to_string(structure)
        view.addModel(pdb_string, "pdb")
    view.setStyle({"cartoon": {"colorscheme": "chain"}})
    for c, i in enumerate(chains_to_align):
        view.setStyle(
            {"model": c}, {"cartoon": {"color": "orange", "opacity": opacity}}
        )
        if isinstance(i, str):
            view.setStyle(
                {"model": c, "chain": i},
                {"cartoon": {"color": "gray", "opacity": opacity}},
            )
        if isinstance(i, list):
            for j in i:
                view.setStyle(
                    {"model": c, "chain": j},
                    {"cartoon": {"color": "gray", "opacity": opacity}},
                )
    view.zoomTo()
    return view


def plot_structure_and_frag_center_of_mass(
    reference_pdb_file,
    fragment_model_pdb_files,
    ref_receptor_chains: str | list | None = None,
    fragment_chain=None,
    opacity=0.5,
):
    # if ref_receptor_chains is None:
    #     ref_receptor_chains = pdb_tools.get_chains_from_structure(reference_pdb_file)[:-1]
    if isinstance(ref_receptor_chains, str):
        ref_receptor_chains = [ref_receptor_chains]
    view = py3Dmol.view(height=600, width=800)
    with open(reference_pdb_file, "r") as f:
        pdb_string = f.read()
    view.addModel(pdb_string, "pdb")
    view.setStyle({"cartoon": {"colorscheme": "chain", "opacity": opacity}})
    if ref_receptor_chains is not None:
        for chain in ref_receptor_chains:
            view.setStyle(
                {"chain": chain}, {"cartoon": {"color": "gray", "opacity": opacity}}
            )
    if fragment_chain is None:
        fragment_chain = pdb_tools.get_chains_from_structure(
            fragment_model_pdb_files[0]
        )[-1]
    for fragment_pdb_file in fragment_model_pdb_files:
        com = pdb_tools.get_chain_center_of_mass(fragment_pdb_file, fragment_chain)
        view.addSphere(
            {
                "center": {"x": com[0], "y": com[1], "z": com[2]},
                "radius": 1,
                "color": "orange",
                "opacity": opacity,
            }
        )
    view.zoomTo()
    return view
