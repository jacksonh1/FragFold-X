import pandas as pd
import fragfold3.config as config
import fragfold3.tools.plotting as plotting
import fragfold3.tools.pdb_tools as pdb_tools
import yaml
import matplotlib.pyplot as plt

plt.style.use("fragfold3.local_data.fragfold3_plotstyle")
import seaborn as sns
from pathlib import Path
import shutil
import fragfold3.tools.chimerax_tools as chimerax_tools


def plot_fragfold_results_wrapper(
    run_name, results_csv_dict, fragfold_guide=config.FRAGFOLD_GUIDE, cols=None, **kwargs
):
    csv_file = results_csv_dict[run_name]
    fragment_source_label = fragfold_guide[run_name]["fragments"]
    receptor_labels = fragfold_guide[run_name]["receptor"]
    df = pd.read_csv(csv_file)
    title = (
        f"fragment source: {fragment_source_label} - receptor(s): {' + '.join(receptor_labels)}"
    )
    fig, ax = plotting.plot_fragfold_results(
        df=df,
        title=title,
        cols=cols,
        **kwargs,
    )
    return fig, ax, df


def plotly_fragfold_results_wrapper(
    run_name,
    results_csv_dict,
    fragfold_guide=config.FRAGFOLD_GUIDE,
    cols=None,
    xcol="fragment_center",
    weight_set=None,
):
    csv_file = results_csv_dict[run_name]
    run_info = fragfold_guide[run_name]
    fragment_source = run_info["fragments"]
    receptors = run_info["receptor"]
    return plotting.plotly_fragfold_results(
        csv_file,
        fragment_source_label=fragment_source,
        receptor_labels=receptors,
        cols=cols,
        xcol=xcol,
        weight_set=weight_set,
    )


def get_pdbs_from_region(
    dfin,
    start,
    end,
    weights="alphafold2_ptm",
    ranks=None,
    root=config.PROJ_ROOT,
    pdb_colname="pdb_file_relative",
):
    if ranks is None:
        ranks = [1, 2, 3, 4, 5]
    if not isinstance(ranks, list):
        ranks = [ranks]
    df = dfin.copy()
    x = df[
        (df["fragment_center"] > start)
        & (df["fragment_center"] < end)
        & (df["weights"] == weights)
        & (df["rank"].isin(ranks))
    ][pdb_colname].to_list()
    return [root / i for i in x]


# def add_aligned_pdb(pdb_file, aligned_foldername="aligned_pdbs", root=config.PROJ_ROOT):
#     pdb_file = root / pdb_file
#     aln_pdb_filename = pdb_file.stem + "-aligned.pdb"
#     aligned_pdb_file = pdb_file.parent / aligned_foldername / aln_pdb_filename
#     if not aligned_pdb_file.exists():
#         raise FileNotFoundError(f"Aligned pdb file {aligned_pdb_file} not found")
#     return aligned_pdb_file.resolve().relative_to(root)


def add_ref_pdb(pdb_file, ref_pdb_name=None, root=config.PROJ_ROOT):
    pdb_file = root / pdb_file
    if ref_pdb_name is None:
        ref_pdb_name = Path(
            config.FRAGFOLD_GUIDE[pdb_file.resolve().parent.parent.parent.name]["reference_pdb"]
        ).name
    ref_pdb_file = pdb_file.parent / ref_pdb_name
    if not ref_pdb_file.exists():
        raise FileNotFoundError(f"Reference pdb file {ref_pdb_file} not found")
    return ref_pdb_file.resolve().relative_to(root)


# def add_aligned_ref_pdb(ref_pdb_file, aligned_foldername="aligned_pdbs", root=config.PROJ_ROOT):
#     ref_pdb_file = root / ref_pdb_file
#     aln_ref_pdb_filename = ref_pdb_file.stem + "-aligned.pdb"
#     aligned_ref_pdb_file = ref_pdb_file.parent / aligned_foldername / aln_ref_pdb_filename
#     if not aligned_ref_pdb_file.exists():
#         raise FileNotFoundError(f"Aligned reference pdb file {aligned_ref_pdb_file} not found")
#     return aligned_ref_pdb_file.resolve().relative_to(root)


def add_fragment_sequence(pdb_file, chain=None, root=config.PROJ_ROOT):
    pdb_file = root / pdb_file
    if chain is None:
        chains = pdb_tools.get_chains_from_structure(pdb_file)
        chain = chains[-1]
    fragment_sequence = pdb_tools.extract_chain_sequence_from_resnames(pdb_file, chain_id=chain)
    return fragment_sequence


def copy_pdbs_to_folder(
    df,
    start,
    end,
    run_name,
    output_folder="./",
    weights="alphafold2_ptm",
    ranks=None,
    root=config.PROJ_ROOT,
    pdb_colname="pdb_file_relative",
):
    pdb_files = get_pdbs_from_region(
        df, start, end, weights=weights, ranks=ranks, root=root, pdb_colname=pdb_colname
    )
    output_folder = Path(output_folder)
    output_folder_path = output_folder / f"{run_name}_{start}_{end}-{weights}"
    output_folder_path.mkdir(parents=True, exist_ok=True)
    # clear output folder
    for f in output_folder_path.glob("*"):
        f.unlink()
    pdb_files.append(root / df[df["weights"] == weights]["aligned_ref_pdb_file_relative"].iloc[0])
    for pdb_file in pdb_files:
        if pdb_file.exists():
            shutil.copy(pdb_file, output_folder_path)
        else:
            print(f"Warning: PDB file {pdb_file} does not exist and will be skipped.")
    print(f"Copied {len(pdb_files)} PDB files to {output_folder_path}")
    return output_folder_path


def create_summary_files(
    df,
    start,
    end,
    run_name,
    export_columns=None,
    output_folder="./",
    weights="alphafold2_ptm",
    ranks=None,
    root=config.PROJ_ROOT,
    pdb_colname="pdb_file_relative",
    fragfold_guide=config.FRAGFOLD_GUIDE,
):
    output_folder_path = copy_pdbs_to_folder(
        df,
        start,
        end,
        run_name,
        output_folder=output_folder,
        weights=weights,
        ranks=ranks,
        root=root,
        pdb_colname=pdb_colname,
    )
    if ranks is None:
        ranks = [1, 2, 3, 4, 5]
    if not isinstance(ranks, list):
        ranks = [ranks]
    if export_columns is None:
        export_columns = [
            i
            for i in df.columns
            if i not in ["aligned_ref_pdb_file_relative", "aligned_pdb_file_relative", "contacts"]
        ]
    fragment_source_label = fragfold_guide[run_name]["fragments"]
    receptor_labels = fragfold_guide[run_name]["receptor"]
    title = (
        f"fragment source: {fragment_source_label} - receptor(s): {' + '.join(receptor_labels)}"
    )
    full_plot_file = output_folder_path / f"{run_name}_{start}_{end}_full.png"
    fig, axes = plotting.plot_fragfold_results(
        df=df,
        title=title,
    )
    axes = plotting.highlight_regions(axes, start, end)
    fig.savefig(full_plot_file, bbox_inches="tight", dpi=300)
    plt.close(fig)
    just_spec_weight_plot_file = output_folder_path / f"{run_name}_{start}_{end}-{weights}.png"
    df_weights = df[df["weights"] == weights].copy()
    fig, axes = plotting.plot_fragfold_results(
        df=df_weights,
        title=title,
    )
    axes = plotting.highlight_regions(axes, start, end)
    fig.savefig(just_spec_weight_plot_file, bbox_inches="tight", dpi=300)
    plt.close(fig)
    x = df[
        (df["fragment_center"] > start)
        & (df["fragment_center"] < end)
        & (df["weights"] == weights)
        & (df["rank"].isin(ranks))
    ].copy()
    x[export_columns].to_csv(
        output_folder_path / f"{output_folder_path.name}-scores.csv", index=False
    )
    chimerax_tools.align_pdbs_chimerax(output_folder_path)
