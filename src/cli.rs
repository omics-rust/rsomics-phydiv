use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use rsomics_phylo_tree::Tree;

use rsomics_phydiv::{Config, Rooted, Weight, run};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-phydiv", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Count table (feature-by-sample TSV); reads stdin when "-" or omitted.
    #[arg(default_value = "-")]
    input: PathBuf,

    /// Newick tree whose tips are the OTU/taxon IDs.
    #[arg(long)]
    tree: PathBuf,

    /// Force the rooted variant (include the root path).
    #[arg(long, conflicts_with = "unrooted")]
    rooted: bool,

    /// Force the unrooted variant (drop branches above the community LCA).
    #[arg(long)]
    unrooted: bool,

    /// Abundance weighting: 0 (unweighted), 1 (fully weighted), or θ in (0,1).
    #[arg(long, default_value = "0")]
    weight: String,

    /// Treat the input table as comma-separated instead of tab-separated.
    #[arg(long, default_value_t = false)]
    csv: bool,

    /// Decimal places in the output.
    #[arg(short = 'p', long, default_value_t = 6)]
    precision: usize,

    /// Output path; writes stdout when "-".
    #[arg(short = 'o', long, default_value = "-")]
    output: String,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let newick = fs::read_to_string(&self.tree)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", self.tree.display())))?;
        let tree = Tree::from_newick(&newick)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", self.tree.display())))?;

        let rooted = if self.rooted {
            Rooted::Rooted
        } else if self.unrooted {
            Rooted::Unrooted
        } else {
            Rooted::Auto
        };

        let cfg = Config {
            delim: if self.csv { ',' } else { '\t' },
            rooted,
            weight: Weight::parse(&self.weight)?,
            precision: self.precision,
        };

        let reader: Box<dyn std::io::BufRead> = if self.input.as_os_str() == "-" {
            Box::new(BufReader::new(std::io::stdin().lock()))
        } else {
            Box::new(BufReader::new(File::open(&self.input).map_err(|e| {
                RsomicsError::InvalidInput(format!("{}: {e}", self.input.display()))
            })?))
        };
        let mut out: Box<dyn Write> = if self.output == "-" && self.common.json {
            Box::new(BufWriter::new(std::io::sink()))
        } else if self.output == "-" {
            Box::new(BufWriter::new(std::io::stdout().lock()))
        } else {
            Box::new(BufWriter::new(
                File::create(&self.output).map_err(RsomicsError::Io)?,
            ))
        };
        run(reader, &mut out, &tree, &cfg)?;
        out.flush().map_err(RsomicsError::Io)
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Generalized phylogenetic alpha-diversity (rooted/unrooted, abundance-weighted PD).",
    origin: Some(Origin {
        upstream: "scikit-bio skbio.diversity.alpha.phydiv",
        upstream_license: "BSD-3-Clause",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.7717/peerj.157"),
    }),
    usage_lines: &[
        "[table.tsv] --tree tree.nwk [--rooted|--unrooted] [--weight 0|1|θ] [-o pd.tsv]",
    ],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "tree",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("path"),
                required: true,
                default: None,
                description: "Newick tree whose tips are the OTU/taxon IDs.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "rooted",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "Force the rooted variant (include the root path).",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "unrooted",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "Force the unrooted variant.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "weight",
                aliases: &[],
                value: Some("<float>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("0"),
                description: "Abundance weighting: 0 unweighted, 1 fully weighted, θ in (0,1).",
                why_default: Some("0 reproduces Faith's PD / uPD."),
            },
            FlagSpec {
                short: None,
                long: "csv",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: Some("false"),
                description: "Parse the table as comma-separated.",
                why_default: None,
            },
            FlagSpec {
                short: Some('p'),
                long: "precision",
                aliases: &[],
                value: Some("<int>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("6"),
                description: "Decimal places in the output.",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("path"),
                required: false,
                default: Some("-"),
                description: "Output path (- for stdout).",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Rooted Faith's PD (auto-detects rooted tree, default weight 0)",
            command: "rsomics-phydiv counts.tsv --tree tree.nwk",
        },
        Example {
            description: "Abundance-weighted unrooted BWPD with θ=0.25",
            command: "rsomics-phydiv counts.tsv --tree tree.nwk --unrooted --weight 0.25",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
