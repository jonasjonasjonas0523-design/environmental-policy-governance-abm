import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({"font.family": "DejaVu Sans", "mathtext.default": "it"})

FIG_W, FIG_H = 10.24, 5.79  # same page size as the supplied Figure1.pdf
fig = plt.figure(figsize=(FIG_W, FIG_H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

EDGE = "#333333"; GRAY = "#e8e8e8"; LIGHT = "#f4f4f4"

def box(x0, y0, x1, y1, text, fs=8.6, fill="white", bold_first=False, lw=1.2, pad=0.35):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                 boxstyle=f"round,pad={pad}", linewidth=lw,
                 edgecolor=EDGE, facecolor=fill, mutation_aspect=FIG_W / FIG_H))
    ax.text((x0 + x1) / 2, (y0 + y1) / 2, text, ha="center", va="center",
            fontsize=fs, linespacing=1.45)

def arrow(p, q, lw=1.3):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=13,
                                 linewidth=lw, color="black", shrinkA=0, shrinkB=0))

def seg(p, q, lw=1.3):
    ax.plot([p[0], q[0]], [p[1], q[1]], color="black", lw=lw, solid_capstyle="round")

# ---- outer governance frame ----
ax.add_patch(FancyBboxPatch((2, 2.5), 96, 94.5, boxstyle="round,pad=0.4",
             linewidth=1.6, edgecolor="#444444", facecolor=LIGHT,
             linestyle=(0, (6, 3)), mutation_aspect=FIG_W / FIG_H))
ax.text(50, 94.0, "Governance environment  (properties of the environment, shared by all instruments)",
        ha="center", va="center", fontsize=10.2, fontweight="bold", color="#222222")
ax.text(50, 90.2, "$p$ detection   $\\cdot$   $c$ consequence credibility   $\\cdot$   "
                  "$k$ implementation leakage   $\\cdot$   $h$ policy horizon",
        ha="center", va="center", fontsize=9.0, color="#444444")

# ---- top row: regulator and realised-consequences boxes ----
box(6, 63, 36, 86,
    "REGULATOR  (one instrument per run)\n"
    "Standard: thresholds $\\bar{E}_i$, penalty $f_S$\n"
    "Charge: rate $\\tau$, evasion penalty $f_T$\n"
    "Voluntary: targets, threat $\\theta_t$", fs=8.8)

box(58, 63, 95, 86,
    "Realised consequences and public signals\n"
    "penalties, charge collection, evaded-charge recovery;\n"
    "verified achievement $z_{i,t}$ updates reputation $R$;\n"
    "recognised shortfall $s_t$ updates threat $\\theta_{t+1}$ (voluntary)", fs=8.4)

# theta feedback: realised -> regulator, routed through the open mid-band
seg((66, 62.4), (66, 58)); seg((66, 58), (21, 58)); arrow((21, 58), (21, 62.4))
ax.text(43.5, 60.1, "$\\theta_{t+1}$ to regulator state (voluntary)",
        ha="center", va="center", fontsize=8.0, style="italic", color="#333333")

# expected-future-consequence feedback: realised -> firms, lower lane of the mid-band
seg((74, 62.4), (74, 54)); seg((74, 54), (30, 54)); arrow((30, 54), (30, 51.8))
ax.text(53.0, 51.4, "expected future consequence  $\\delta h\\,V_{i,t+1}$;   reputation $R$",
        ha="center", va="center", fontsize=8.0, style="italic", color="#333333")

# regulator -> firms
arrow((16, 62.2), (16, 51.6))
ax.text(14.4, 56.8, "policy rules and\ninstrument state",
        ha="right", va="center", fontsize=7.8, style="italic", color="#333333")

# ---- middle row: firms and the enforcement chain ----
box(6, 28, 36, 51,
    "N HETEROGENEOUS FIRMS\n"
    "abatement cost $C_i = \\frac{1}{2}\\gamma_i a^2$,  baseline $E_i^0$\n"
    "minimise  $J = C + L + Q - R + \\delta h\\,V$\n"
    "choose: abatement $a_{i,t}$, report $\\tilde{E}_{i,t}$,\n"
    "participation $d_{i,t}$", fs=8.6)

box(40, 32, 54, 47, "True emissions $E_{i,t}$\nand reports $\\tilde{E}_{i,t}$", fs=8.6)
box(58, 32, 70, 47, "Monitoring /\nverification\n(prob. $p$)", fs=8.6, fill=GRAY)
box(73.5, 32, 84.5, 47, "Formal\nimplementation\n(prob. $c$)", fs=8.6, fill=GRAY)
box(88, 32, 96, 47, "Leakage-\nadjusted\nrealisation\n$(1-k)$", fs=8.4, fill=GRAY)

arrow((36.8, 39.5), (39.2, 39.5))
arrow((54.8, 39.5), (57.2, 39.5))
arrow((70.8, 39.5), (72.7, 39.5))
arrow((85.3, 39.5), (87.2, 39.5))
# realised consequences fed back up from the leakage-adjusted stage
arrow((92, 47.9), (92, 62.2))

# ---- comparison-design note band (previously the overlapping element) ----
box(6, 13.5, 94, 22.5,
    "Same generated firm population;  shared governance functions;  "
    "matched stochastic streams where structurally comparable.\n"
    "Principal estimand:  paired difference  "
    "$\\Delta_{jl}(g) = Y_j(g) - Y_l(g)$  at each  $g = (p,\\,c,\\,k,\\,h)$",
    fs=8.4, fill="white")

# ---- outcome logs strip ----
box(14, 4.2, 86, 10.2,
    "Outcome logs:   abatement $A$;   social resource cost $C_{social}$ "
    "$(\\kappa_M,\\,\\kappa_A,\\,\\kappa_I)$;   fiscal transfers $T_{fiscal}$  (separate)",
    fs=8.4, fill="white")

fig.savefig("Figure1_fixed.pdf", bbox_inches=None)
fig.savefig("figs/Figure1.png", dpi=300, bbox_inches=None)
print("figure written")
