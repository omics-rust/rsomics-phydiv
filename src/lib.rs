use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Write};

use rsomics_common::{Result, RsomicsError};
use rsomics_phylo_tree::{NodeId, Tree};

pub struct CountTable {
    pub feature_ids: Vec<String>,
    pub sample_names: Vec<String>,
    /// Column-major: one count vector per sample, indexed by feature row.
    pub columns: Vec<Vec<u64>>,
}

impl CountTable {
    /// # Errors
    /// Errors on a missing header, a ragged row, or a non-integer count.
    pub fn parse<R: BufRead>(reader: R, delim: char) -> Result<CountTable> {
        let mut lines = reader.lines();
        let header = loop {
            match lines.next() {
                Some(line) => {
                    let line = line.map_err(RsomicsError::Io)?;
                    if line.trim().is_empty() || line.starts_with('#') {
                        continue;
                    }
                    break line;
                }
                None => return Err(RsomicsError::InvalidInput("empty count table".into())),
            }
        };
        let sample_names: Vec<String> = header
            .split(delim)
            .skip(1)
            .map(|s| s.trim().to_string())
            .collect();
        if sample_names.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "header has no sample columns (need feature-ID column + ≥1 sample)".into(),
            ));
        }
        let n = sample_names.len();
        let mut feature_ids = Vec::new();
        let mut seen_features = HashSet::new();
        let mut columns: Vec<Vec<u64>> = vec![Vec::new(); n];
        for (row_idx, line) in lines.enumerate() {
            let line = line.map_err(RsomicsError::Io)?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let mut fields = line.split(delim);
            let feature = fields.next().unwrap_or("").trim().to_string();
            if !seen_features.insert(feature.clone()) {
                return Err(RsomicsError::InvalidInput(format!(
                    "duplicate taxon '{feature}' in the count table; all taxa must be unique"
                )));
            }
            let mut seen = 0usize;
            for (col, field) in fields.enumerate() {
                if col >= n {
                    return Err(RsomicsError::InvalidInput(format!(
                        "row {} (feature '{feature}') has more columns than the header",
                        row_idx + 2
                    )));
                }
                let count: u64 = field.trim().parse().map_err(|_| {
                    RsomicsError::InvalidInput(format!(
                        "row {} (feature '{feature}'), sample '{}': '{}' is not a non-negative integer count",
                        row_idx + 2,
                        sample_names[col],
                        field.trim()
                    ))
                })?;
                columns[col].push(count);
                seen += 1;
            }
            if seen != n {
                return Err(RsomicsError::InvalidInput(format!(
                    "row {} (feature '{feature}') has {seen} count columns, header has {n}",
                    row_idx + 2
                )));
            }
            feature_ids.push(feature);
        }
        Ok(CountTable {
            feature_ids,
            sample_names,
            columns,
        })
    }
}

/// Whether the root branch (and the path above the community LCA) counts.
/// `Auto` follows scikit-bio: rooted iff the root is bifurcating.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Rooted {
    Auto,
    Rooted,
    Unrooted,
}

/// Abundance weighting of branch contributions — the BWPD family of McCoy & Matsen 2013.
/// `Theta` is the partial-weighting exponent θ ∈ (0, 1); θ=0 is unweighted, θ=1 is `Full`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Weight {
    Unweighted,
    Full,
    Theta(f64),
}

impl Weight {
    /// Parse `--weight`: `0` → unweighted, `1` → full, a float in (0,1) → θ.
    ///
    /// # Errors
    /// Errors when the value is not a number in [0, 1].
    pub fn parse(s: &str) -> Result<Weight> {
        let v: f64 = s
            .trim()
            .parse()
            .map_err(|_| RsomicsError::InvalidInput(format!("--weight '{s}' is not a number")))?;
        if !(0.0..=1.0).contains(&v) {
            return Err(RsomicsError::InvalidInput(
                "--weight must be within [0, 1]".into(),
            ));
        }
        Ok(if v == 0.0 {
            Weight::Unweighted
        } else if v == 1.0 {
            Weight::Full
        } else {
            Weight::Theta(v)
        })
    }
}

pub struct Config {
    pub delim: char,
    pub rooted: Rooted,
    pub weight: Weight,
    pub precision: usize,
}

/// Tree resolved for phydiv: branch length per node (every non-root node must
/// carry one; the root's is 0.0 since no branch sits above it), a tip-name →
/// node map, the node ids in postorder, and whether the root is bifurcating.
struct PhyTree {
    branch_length: Vec<f64>,
    children: Vec<Vec<NodeId>>,
    tip_index: HashMap<String, NodeId>,
    postorder: Vec<NodeId>,
    n_nodes: usize,
    root_bifurcating: bool,
}

impl PhyTree {
    fn build(tree: &Tree) -> Result<PhyTree> {
        let n_nodes = tree.nodes.len();
        let mut branch_length = vec![0.0f64; n_nodes];
        let mut children = vec![Vec::new(); n_nodes];
        let mut tip_index = HashMap::new();
        for node in &tree.nodes {
            children[node.id] = node.children.clone();
            match node.branch_length {
                Some(bl) => branch_length[node.id] = bl,
                None if node.id != tree.root => {
                    return Err(RsomicsError::InvalidInput(
                        "all non-root nodes in the tree must have a branch length".into(),
                    ));
                }
                None => {}
            }
            if node.children.is_empty() && node.id != tree.root {
                let name = node
                    .name
                    .as_deref()
                    .ok_or_else(|| RsomicsError::InvalidInput("a tip has no name".into()))?;
                if tip_index.insert(name.to_string(), node.id).is_some() {
                    return Err(RsomicsError::InvalidInput(format!(
                        "duplicate tip name '{name}' in the tree"
                    )));
                }
            }
        }

        let mut postorder = Vec::with_capacity(n_nodes);
        let mut stack = vec![(tree.root, false)];
        while let Some((id, visited)) = stack.pop() {
            if visited {
                postorder.push(id);
            } else {
                stack.push((id, true));
                for &c in &children[id] {
                    stack.push((c, false));
                }
            }
        }

        Ok(PhyTree {
            branch_length,
            children,
            tip_index,
            postorder,
            n_nodes,
            root_bifurcating: tree.nodes[tree.root].children.len() == 2,
        })
    }

    fn is_rooted(&self, rooted: Rooted) -> bool {
        match rooted {
            Rooted::Rooted => true,
            Rooted::Unrooted => false,
            Rooted::Auto => self.root_bifurcating,
        }
    }

    /// Fill `cbn[id]` with the total abundance of taxa descending from each node,
    /// accumulated up the tree in postorder.
    fn accumulate(&self, tip_counts: &[(NodeId, u64)], cbn: &mut [f64]) {
        cbn.iter_mut().for_each(|c| *c = 0.0);
        for &(tip, c) in tip_counts {
            cbn[tip] = c as f64;
        }
        for &id in &self.postorder {
            if !self.children[id].is_empty() {
                cbn[id] = self.children[id].iter().map(|&c| cbn[c]).sum();
            }
        }
    }
}

/// Generalized phylogenetic diversity from per-node descendant abundances.
/// McCoy & Matsen 2013 BWPD: unrooted drops branches above the LCA and folds each
/// branch to its abundance balance `2·min(p, 1−p)`; θ tempers the weighting.
fn diversity(pt: &PhyTree, cbn: &[f64], rooted: bool, weight: Weight) -> f64 {
    let total = cbn.iter().copied().fold(0.0f64, f64::max);
    if total == 0.0 {
        return 0.0;
    }
    let mut sum = 0.0;
    for (id, &c) in cbn.iter().enumerate() {
        let factor = match weight {
            Weight::Unweighted => {
                if c > 0.0 && (rooted || c < total) {
                    1.0
                } else {
                    0.0
                }
            }
            _ => {
                let mut frac = c / total;
                if !rooted {
                    frac = 2.0 * frac.min(1.0 - frac);
                }
                match weight {
                    Weight::Theta(theta) => frac.powf(theta),
                    _ => frac,
                }
            }
        };
        sum += pt.branch_length[id] * factor;
    }
    sum
}

pub fn run<R: BufRead, W: Write>(reader: R, out: &mut W, tree: &Tree, cfg: &Config) -> Result<()> {
    let table = CountTable::parse(reader, cfg.delim)?;
    let pt = PhyTree::build(tree)?;
    let rooted = pt.is_rooted(cfg.rooted);

    let row_tip: Vec<NodeId> = table
        .feature_ids
        .iter()
        .map(|taxon| {
            pt.tip_index.get(taxon).copied().ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "taxon '{taxon}' from the count table is not a tip in the tree"
                ))
            })
        })
        .collect::<Result<_>>()?;

    writeln!(out, "sample\tphydiv").map_err(RsomicsError::Io)?;
    let mut tip_counts = Vec::new();
    let mut cbn = vec![0.0f64; pt.n_nodes];
    for (col, sample) in table.sample_names.iter().enumerate() {
        tip_counts.clear();
        for (row, &c) in table.columns[col].iter().enumerate() {
            if c > 0 {
                tip_counts.push((row_tip[row], c));
            }
        }
        pt.accumulate(&tip_counts, &mut cbn);
        let value = diversity(&pt, &cbn, rooted, cfg.weight);
        writeln!(out, "{sample}\t{value:.*}", cfg.precision).map_err(RsomicsError::Io)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn doc_tree() -> Tree {
        Tree::from_newick("((a:1,b:2)c:0.5,(d:1,e:1)f:1)root;").unwrap()
    }

    fn phy(tree: &Tree, table: &str, rooted: Rooted, weight: Weight) -> f64 {
        let cfg = Config {
            delim: '\t',
            rooted,
            weight,
            precision: 12,
        };
        let mut out = Vec::new();
        run(std::io::Cursor::new(table), &mut out, tree, &cfg).unwrap();
        String::from_utf8(out)
            .unwrap()
            .lines()
            .nth(1)
            .unwrap()
            .split_once('\t')
            .unwrap()
            .1
            .parse()
            .unwrap()
    }

    const T1: &str = "feature\tu\na\t1\nb\t0\nd\t3\ne\t2\n";

    #[test]
    fn unweighted_matches_faith() {
        assert!((phy(&doc_tree(), T1, Rooted::Rooted, Weight::Unweighted) - 4.5).abs() < 1e-12);
    }

    #[test]
    fn rooted_full_weight() {
        let v = phy(&doc_tree(), T1, Rooted::Rooted, Weight::Full);
        assert!((v - 1.916_666_666_666_666_5).abs() < 1e-12);
    }

    #[test]
    fn unrooted_full_weight() {
        let v = phy(&doc_tree(), T1, Rooted::Unrooted, Weight::Full);
        assert!((v - 2.5).abs() < 1e-12);
    }

    #[test]
    fn rooted_theta_quarter() {
        let v = phy(&doc_tree(), T1, Rooted::Rooted, Weight::Theta(0.25));
        assert!((v - 3.514_589_549_479_082).abs() < 1e-12);
    }

    #[test]
    fn unrooted_theta_half() {
        let v = phy(&doc_tree(), T1, Rooted::Unrooted, Weight::Theta(0.5));
        assert!((v - 3.259_872_253_901_790_4).abs() < 1e-12);
    }

    #[test]
    fn auto_bifurcating_is_rooted() {
        let auto = phy(&doc_tree(), T1, Rooted::Auto, Weight::Full);
        let rooted = phy(&doc_tree(), T1, Rooted::Rooted, Weight::Full);
        assert_eq!(auto, rooted);
    }

    #[test]
    fn auto_trifurcating_is_unrooted() {
        let tree = Tree::from_newick("(a:1,b:2,c:3)root;").unwrap();
        let table = "feature\ts\na\t2\nb\t3\nc\t0\n";
        let auto = phy(&tree, table, Rooted::Auto, Weight::Full);
        let unrooted = phy(&tree, table, Rooted::Unrooted, Weight::Full);
        assert_eq!(auto, unrooted);
    }

    #[test]
    fn rooted_vs_unrooted_differ_on_subset() {
        let tree = Tree::from_newick("(((a:1,b:2)g:3,c:1.5)h:0.7,(d:1,e:1)f:1)root;").unwrap();
        let table = "feature\ts\na\t5\nb\t4\nc\t0\nd\t0\ne\t0\n";
        let r = phy(&tree, table, Rooted::Rooted, Weight::Unweighted);
        let u = phy(&tree, table, Rooted::Unrooted, Weight::Unweighted);
        assert!((r - 6.7).abs() < 1e-12);
        assert!((u - 3.0).abs() < 1e-12);
    }

    #[test]
    fn empty_sample_is_zero() {
        let table = "feature\tz\na\t0\nb\t0\nd\t0\ne\t0\n";
        assert_eq!(phy(&doc_tree(), table, Rooted::Rooted, Weight::Full), 0.0);
    }

    #[test]
    fn weight_parses() {
        assert_eq!(Weight::parse("0").unwrap(), Weight::Unweighted);
        assert_eq!(Weight::parse("1").unwrap(), Weight::Full);
        assert_eq!(Weight::parse("0.25").unwrap(), Weight::Theta(0.25));
        assert!(Weight::parse("1.5").is_err());
        assert!(Weight::parse("x").is_err());
    }

    fn err_cfg() -> Config {
        Config {
            delim: '\t',
            rooted: Rooted::Auto,
            weight: Weight::Unweighted,
            precision: 6,
        }
    }

    fn run_err(tree: &Tree, table: &str) -> bool {
        let mut out = Vec::new();
        run(std::io::Cursor::new(table), &mut out, tree, &err_cfg()).is_err()
    }

    #[test]
    fn unknown_taxon_rejected() {
        assert!(run_err(&doc_tree(), "feature\tx\na\t1\nzzz\t1\n"));
    }

    #[test]
    fn missing_branch_length_rejected() {
        let tree = Tree::from_newick("((a:1,b)c:0.5,(d:1,e:1)f:1)root;").unwrap();
        assert!(run_err(&tree, "feature\tx\na\t1\nb\t1\n"));
    }

    #[test]
    fn duplicate_taxon_rejected() {
        assert!(run_err(&doc_tree(), "feature\tx\na\t1\na\t2\n"));
    }

    #[test]
    fn single_node_tree_has_no_tip() {
        let tree = Tree::from_newick("a;").unwrap();
        assert!(run_err(&tree, "feature\tx\na\t1\n"));
    }
}
