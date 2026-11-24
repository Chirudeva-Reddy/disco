import streamlit as st

from app.utils import constants
from app.utils import helper



#region Config
st.set_page_config(
    page_title=f"Help | {constants.APP_NAME}",
    initial_sidebar_state="expanded",
    # layout="centered",
)

st.html(f"""
<h1 align="center" style="font-size: 48px;">
    Help |
    <img src="{constants.ASSET_LOGO_DATAURI}" alt="" width="64" height="64"/>
    {constants.APP_NAME}
</h1>
""")
st.caption("A Consolidated database of Disorder, Patterns and Functional annotations of Transcription Factors", text_alignment="center")

#endregion

st.divider()

st.sidebar.write(f"""
## **[Definitions](#definitions)**
\n:material/arrow_right_alt: [TFClass hierarchy](#tfclass_hierarchy)
\n:material/arrow_right_alt: [Vagueness metric, Vagueness penalty](#vagueness)
\n:material/arrow_right_alt: [Shuffle rank](#shuffle_rank)
\n:material/arrow_right_alt: [Shuffle K comparison](#shuffle_k)

## **Pages**
\n:material/home: [Home Page](#home_page)
\n:material/view_list: [TF Browser](#tf_browser)
\n:material/visibility: [TF Viewer](#tf_viewer)
\n:material/regular_expression: [Pattern Explorer](#pattern_explorer)
""")



#region Intro

st.header(":primary[:material/dictionary:] **Definitions**", anchor="definitions")

helper.container_indent().markdown(f"""
- :primary[**TF**]: Transcription Factor.
- :primary[**DBD** / **AD**]: DNA Binding Domain / Activation Domain. Note that
  since the annotation of DBD regions can never be exhaustive, a general
  consensus is assumed to consider non-DBD regions as Activation Domains (AD).
- :primary[**ELM Pattern**]: Eukaryotic Linear Motif Pattern/Motif. These are
  short, conserved sequences within proteins that are recognized by specific
  protein domains. These patterns are denoted as Regular Expressions.
- :primary[**RegEx**]: Regular Expression. A sequence of characters that defines
  a search pattern, used for pattern-matching within strings.
- File formats used for data representation.
  - :primary[**CSV**]: Comma Separated Values (tabular format)
  - :primary[**TSV**]: Tab Separated Values (tabular format)
  - :primary[**HTML**]: HyperText Markup Language (web page format)
  - :primary[**FASTA**]: Plaintext-based format for representing protein sequences
""")

st.subheader("TFClass hierarchy", anchor="tfclass_hierarchy")

helper.container_indent().markdown(f"""
The TFClass resource classifies transcription factors into a hierarchical
structure based on their DNA-binding domains and other functional
characteristics. There are 5 levels of classification
([source](http://www.edgar-wingender.de/TFClass_schema.html)):

| Level | Rank denomination | Definition                           | Example                                                       |
| ----: | ----------------- | ------------------------------------ | ------------------------------------------------------------- |
|     1 | Superclass        | General topology of the DBD          | Zinc-coordinating DBDs (Superclass 2)                         |
|     2 | Class             | Structural blueprint of the DBD      | Nuclear receptors with C4 zinc fingers (Class 2.1)            |
|     3 | Family            | Sequence and functional similarities | Thyroid hormone receptor-related factors (NR1) (Family 2.1.2) |
|     4 | Subfamily         | Sequence-based subgroupings          | Retinoic acid receptors (NR1B) (Subfamily 2.1.2.1)            |
|     5 | Genus             | TF gene                              | RAR-α (Genus 2.1.2.1.1)                                       |

For more details, refer to the [TFClass
resource](http://tfclass.bioinf.med.uni-goettingen.de/index.jsf).
""")

st.subheader("Vagueness metric, Vagueness penalty", anchor="vagueness")

helper.container_indent().markdown("""
:primary[**Vagueness**] is a novel metric introduced in this tool, which may be
used to rank patterns/motifs based on the number of substitutions and wildcards
that the pattern allows. This allows users to :primary[prioritize simpler
patterns] with higher matches, allowing to focus on motifs which are
:primary[potentially biologically relevant], when analyzing TFs.

Patterns with a higher vagueness penalty are considered more specific, while
those with a lower vagueness penalty are considered more general or vague.\\
For example, the [pattern](http://elm.eu.org/elms/ELME000450)
`..FF[^P]{0,2}[KR]{1,2}[^P]{0,4}` is quite vague, since it contains many
wildcards (`.`), and allows for an already-vague character class (`[^P]`) to
have a variable number of repeats (`{0,4}`). Thus, the vagueness score is high
for this pattern.\\
In contrast, the [pattern](http://elm.eu.org/elms/ELME000170) `[FY][DEP]WM` is
less vague, since it contains fixed amino acids (`W`, `M`) and the character
classes (`[FY]`, `[DEP]`) have very less variability, allowing for 2-3
substitutions. Thus, the vagueness score is low for this pattern.

:primary[**Vagueness penalty**] is a subjective, user-selected penalty applied
to the parts of the regex pattern which are highly variable in nature, such as
wildcards and huge character classes. A value of 1.0 indicates no penalty, and
thus the user can set a value greater than 1.0 to introduce the penalty.

It is to be noted that this is only used to **rank** the patterns for
convenience purposes, and it doesn't affect the actual pattern matching, and it
doesn't discard any patterns. If desired, the user can still click on a column,
such as "Observed matches" to sort the patterns based on the number of matches
instead.
""")

st.subheader("Shuffle rank", anchor="shuffle_rank")

helper.container_indent().markdown("""
:primary[**Shuffle rank**] orders the patterns within one TF by how rarely a
randomly reordered copy of that TF matches them as often as the real sequence.

1. Residue 1 is kept in place, and residues 2 onwards are randomly reordered
   ([Wregex 3.0](https://pmc.ncbi.nlm.nih.gov/articles/PMC11318550/)).
2. The non-overlapping matches of each pattern are counted in the reordered copy.
3. This is repeated $K = 200{,}000$ times, and the :primary[**Shuffle p**] is
   $(1 + b) / (K + 1)$, where $b$ is the number of copies with at least as many
   matches as the real sequence
   ([Phipson and Smyth, 2010](https://doi.org/10.2202/1544-6115.1585);
   [Hemerik and Goeman, 2018](https://doi.org/10.1007/s11749-017-0571-1)).
4. The smallest Shuffle p gets rank 1. Equal values share a rank.

Repeated shuffles of a sequence to count ELM matches are also used by
[mimicINT](https://pmc.ncbi.nlm.nih.gov/articles/PMC11982803/) and
[Hagai et al., 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4089993/).
Using a whole-protein shuffle to order the patterns within a TF is our
application of these sources.

It is to be noted that a high rank means that the count is unusual when the
residues of this TF are reordered. It doesn't show that a match is functional.
The values are not adjusted for multiple testing, and the Vagueness penalty
doesn't affect them.
""")

st.subheader("Shuffle K comparison", anchor="shuffle_k")

helper.container_indent().markdown("""
The :primary[**Shuffle K Comparison**] page shows how the shuffle rank of each
TF changes with the number of shuffles $K$. One stream of $200{,}000$ shuffles
per TF is read off at $K = 10{,}000$, $50{,}000$, $100{,}000$ and $200{,}000$, so
the smaller $K$ are the first shuffles of the larger one. Every pair of $K$ is
compared, and the page opens with a plain-language guide to each measure.

- The :primary[**top 10%**] of a TF are the classes whose rank is at most
  $\\lceil 0.1 n \\rceil$, where $n$ is the number of classes in the TF. Classes
  tied at the cut are all kept.
- :primary[**Top-10% overlap**] is $|A \\cap B| / |A \\cup B|$ between the top
  sets at two $K$ ([Jaccard, 1912](https://doi.org/10.1111/j.1469-8137.1912.tb05611.x)).
- :primary[**Spearman ρ**] compares the full rankings (1 = same order)
  ([Spearman, 1904](https://doi.org/10.2307/1412159)).
- :primary[**At floor**] counts classes which no shuffle reached. They all have
  the smallest possible $p = 1/(K+1)$ and so tie; a larger $K$ can separate them.

Comparing $K$ this way is our application of the sources in
[Shuffle rank](#shuffle_rank).
""")

#endregion

st.divider()

#region Pages
st.header(":primary[:material/pageview:] **Pages**", anchor="pages")

st.markdown(f"""
{constants.APP_NAME} is organized into a few pages, and the user can perform the
following workflow:
1. Browse and discover available TFs
2. Analyze a particular TF in detail
3. Explore patterns/motifs across a selection of TFs or the entire dataset
""")

st.info("""Throughout the app, there are :primary[:material/help:] help icons
near inputs and widgets, which give a brief explanation on how to use them.""",
icon=":material/info:")

st.markdown("Below is an overview of each page and its functionality:")

st.subheader(":material/home: Home Page", anchor="home_page")
helper.container_indent().markdown(
f"""
The user can start search for a TF in the search-box using [UniProt
accession](https://www.uniprot.org/help/accession_numbers), Genus number, or
Genus name (derived from the [TFClass
Resource](http://tfclass.bioinf.med.uni-goettingen.de/index.jsf)), and use the
:primary[:material/search: Go To Viewer] button to explore the details in the
:primary[:material/visibility: TF Viewer].

Or click on one of the :primary[example TFs] to explore how
{constants.APP_NAME} works.\\
Or go to :primary[:material/view_list: TF Browser] to see all available TFs.
""")

st.subheader(":material/view_list: TF Browser", anchor="tf_browser")
helper.container_indent().markdown(constants.CONTENT_HELP_TF_BROWSER)

st.subheader(":material/visibility: TF Viewer", anchor="tf_viewer")
helper.container_indent().markdown(constants.CONTENT_HELP_TF_VIEWER)

st.subheader(":material/regular_expression: Pattern Explorer", anchor="pattern_explorer")
helper.container_indent().markdown(constants.CONTENT_HELP_PATTERN_EXPLORER(is_on_help_page=True))

#endregion

st.divider()

#region Contact
st.header(":primary[:material/contact_mail:] **Contact**", anchor="contact")
st.markdown(constants.CONTENT_CONTACT)

#endregion
