import { useState } from "react";

// ─────────────────────────────────────────────────────────────────────────────
// MODEL PARAMETERS (post-estimation shrinkage)
// ─────────────────────────────────────────────────────────────────────────────

// CVD Model — Supplemental Table 3
const CVD = {
  alpha:         -6.046429053,
  b_leeftijd:     0.040672727,
  b_vrouw:       -0.234111177,
  b_dmduur:       0.013062752,
  b_sbp:          0.005814221,
  b_ldl:          0.082287009,    // mmol/L
  b_hba1c:        0.012209026,    // mmol/mol
  b_micro:        0.437359313,
  b_macro:        0.738916137,
  b_egfr_lt40:   -0.404318528,    // log2(eGFR), leeftijd < 40
  b_egfr_ge40:   -0.345596046,    // log2(eGFR), leeftijd ≥ 40
  b_roken:        0.204224209,
  b_geenBeweging: 0.229279688,
};

// ESKD Core Model — Supplemental Table 4
const ESKD = {
  alpha:    4.783935,
  b_leeftijd: -0.015820,
  b_man:       0.229484,
  b_dmduur:    0.009617,
  b_egfr:     -1.975550,
  b_micro:     0.063436,
  b_macro:     0.482441,
  b_sbp:       0.006247,
  b_hba1c:     0.006987,
  b_roken:     0.265174,
  b_cvd:       0.243471,
  // Sterfte (niet-ESKD)
  sigma:              -10.193500,
  g_man:                1.746154,
  g_dmduur:             0.008748,
  g_egfr:              -0.044020,
  g_micro:              0.380168,
  g_macro:              0.666386,
  g_sbp:               -0.011110,
  g_hba1c:              0.007618,
  g_roken:              0.679023,
  g_cvd:                0.601755,
  g_leeftijd_vrouw:     0.091864,
  g_leeftijd_man:       0.069326,
};

// ─────────────────────────────────────────────────────────────────────────────
// ACR (mg/mmol) → albuminurie categorie
// Grenswaarden gebaseerd op Nederlandse NfN/NIV richtlijn:
//   Normoalbuminurie : vrouw <3,5 mg/mmol  |  man <2,5 mg/mmol
//   Micro-albuminurie: vrouw 3,5–35 mg/mmol |  man 2,5–25 mg/mmol
//   Macro-albuminurie: vrouw ≥35 mg/mmol   |  man ≥25 mg/mmol
// ─────────────────────────────────────────────────────────────────────────────
function acrNaarCategorie(acr, vrouw) {
  const drempelNormo = vrouw ? 3.5 : 2.5;
  const drempelMicro = vrouw ? 35.0 : 25.0;
  if (acr < drempelNormo) return "normo";
  if (acr < drempelMicro) return "micro";
  return "macro";
}

function albLabel(acr, vrouw) {
  const cat = acrNaarCategorie(acr, vrouw);
  if (cat === "normo") return { tekst: "Normaal", kleur: "#2e7d32" };
  if (cat === "micro") return { tekst: "Micro-albuminurie", kleur: "#e65100" };
  return { tekst: "Macro-albuminurie", kleur: "#b71c1c" };
}

// ─────────────────────────────────────────────────────────────────────────────
// BEREKENINGSFUNCTIES
// ─────────────────────────────────────────────────────────────────────────────
function berekenCVDRisico(inv, jaren) {
  const alb = acrNaarCategorie(inv.acr, inv.vrouw);
  const micro = alb === "micro" ? 1 : 0;
  const macro = alb === "macro" ? 1 : 0;
  const egfrTerm = inv.leeftijd < 40
    ? CVD.b_egfr_lt40 * Math.log2(inv.egfr)
    : CVD.b_egfr_ge40 * Math.log2(inv.egfr);

  const logRate = CVD.alpha
    + CVD.b_leeftijd     * inv.leeftijd
    + CVD.b_vrouw        * (inv.vrouw ? 1 : 0)
    + CVD.b_dmduur       * inv.dmDuur
    + CVD.b_sbp          * inv.sbp
    + CVD.b_ldl          * inv.ldl
    + CVD.b_hba1c        * inv.hba1c
    + CVD.b_micro        * micro
    + CVD.b_macro        * macro
    + egfrTerm
    + CVD.b_roken        * (inv.roken ? 1 : 0)
    + CVD.b_geenBeweging * (inv.bewegen ? 0 : 1);

  return (1 - Math.exp(-Math.exp(logRate) * jaren)) * 100;
}

function berekenESKDRisico(inv, jaren) {
  const alb = acrNaarCategorie(inv.acr, inv.vrouw);
  const micro = alb === "micro" ? 1 : 0;
  const macro = alb === "macro" ? 1 : 0;
  const man = inv.vrouw ? 0 : 1;

  const lE = Math.exp(
    ESKD.alpha
    + ESKD.b_leeftijd * inv.leeftijd
    + ESKD.b_man      * man
    + ESKD.b_dmduur   * inv.dmDuur
    + ESKD.b_egfr     * Math.log2(inv.egfr)
    + ESKD.b_micro    * micro
    + ESKD.b_macro    * macro
    + ESKD.b_sbp      * inv.sbp
    + ESKD.b_hba1c    * inv.hba1c
    + ESKD.b_roken    * (inv.roken ? 1 : 0)
    + ESKD.b_cvd      * (inv.eerderCVD ? 1 : 0)
  );

  const lD = Math.exp(
    ESKD.sigma
    + ESKD.g_man      * man
    + ESKD.g_dmduur   * inv.dmDuur
    + ESKD.g_egfr     * Math.log2(inv.egfr)
    + ESKD.g_micro    * micro
    + ESKD.g_macro    * macro
    + ESKD.g_sbp      * inv.sbp
    + ESKD.g_hba1c    * inv.hba1c
    + ESKD.g_roken    * (inv.roken ? 1 : 0)
    + ESKD.g_cvd      * (inv.eerderCVD ? 1 : 0)
    + (inv.vrouw ? ESKD.g_leeftijd_vrouw : ESKD.g_leeftijd_man) * inv.leeftijd
  );

  const tot = lE + lD;
  return (lE / tot) * (1 - Math.exp(-tot * jaren)) * 100;
}

// ─────────────────────────────────────────────────────────────────────────────
// RISICO KLEUR
// ─────────────────────────────────────────────────────────────────────────────
function risicoKleur(pct) {
  if (pct < 10) return { bg: "#e8f5e9", text: "#2e7d32", label: "Laag" };
  if (pct < 20) return { bg: "#fff8e1", text: "#f57f17", label: "Matig" };
  if (pct < 40) return { bg: "#fff3e0", text: "#e65100", label: "Hoog" };
  return { bg: "#fce4ec", text: "#b71c1c", label: "Zeer hoog" };
}

// ─────────────────────────────────────────────────────────────────────────────
// METER COMPONENT
// ─────────────────────────────────────────────────────────────────────────────
function Meter({ pct, label, sublabel }) {
  const begrensd = Math.min(Math.max(pct, 0), 100);
  const kleur = risicoKleur(pct);
  const omtrek = Math.PI * 54;
  const streep = (begrensd / 100) * omtrek;

  return (
    <div style={{
      display: "flex", flexDirection: "column", alignItems: "center",
      background: kleur.bg, borderRadius: 16, padding: "18px 14px 12px",
      minWidth: 118, flex: 1,
      border: `1.5px solid ${kleur.text}28`,
      transition: "background 0.4s",
    }}>
      <svg width={124} height={68} viewBox="0 0 124 68">
        <path d="M 7 63 A 55 55 0 0 1 117 63"
          fill="none" stroke="#e0e0e0" strokeWidth={9} strokeLinecap="round" />
        <path d="M 7 63 A 55 55 0 0 1 117 63"
          fill="none" stroke={kleur.text} strokeWidth={9} strokeLinecap="round"
          strokeDasharray={`${streep} ${omtrek}`}
          style={{ transition: "stroke-dasharray 0.9s cubic-bezier(0.4,0,0.2,1)" }} />
        <text x="62" y="56" textAnchor="middle"
          style={{ fontSize: 20, fontWeight: 800, fill: kleur.text, fontFamily: "'DM Serif Display', serif" }}>
          {pct.toFixed(1)}%
        </text>
      </svg>
      <div style={{ fontSize: 11, fontWeight: 700, color: kleur.text, letterSpacing: "0.07em", textTransform: "uppercase", marginTop: 2 }}>
        {kleur.label}
      </div>
      <div style={{ fontSize: 12, fontWeight: 600, color: "#2d3748", marginTop: 3, textAlign: "center" }}>{label}</div>
      <div style={{ fontSize: 11, color: "#a0aec0", marginTop: 1 }}>{sublabel}</div>
    </div>
  );
}

// ─────────────────────────────────────────────────────────────────────────────
// FORMULIER HULPCOMPONENTEN
// ─────────────────────────────────────────────────────────────────────────────
function Veld({ label, hint, fout, children }) {
  return (
    <div style={{ marginBottom: 12 }}>
      <label style={{
        display: "block", fontSize: 11, fontWeight: 700, color: "#4a5568",
        letterSpacing: "0.07em", textTransform: "uppercase", marginBottom: 5,
      }}>
        {label}
        {hint && <span style={{ fontWeight: 400, textTransform: "none", color: "#a0aec0", marginLeft: 5 }}>{hint}</span>}
      </label>
      {children}
      {fout && <div style={{ color: "#e53935", fontSize: 11, marginTop: 3 }}>⚠ {fout}</div>}
    </div>
  );
}

const invoerStijl = {
  width: "100%", padding: "9px 11px", borderRadius: 8,
  border: "1.5px solid #dde1e7", fontSize: 14, color: "#1a202c",
  background: "#f9fafb", outline: "none", boxSizing: "border-box",
  fontFamily: "inherit", transition: "border-color 0.2s",
};

function Schakelaar({ waarde, onChange, opties }) {
  return (
    <div style={{ display: "flex", borderRadius: 8, overflow: "hidden", border: "1.5px solid #dde1e7" }}>
      {opties.map(opt => (
        <button key={opt.waarde} onClick={() => onChange(opt.waarde)} style={{
          flex: 1, padding: "9px 4px", fontSize: 13, fontWeight: 600,
          cursor: "pointer", border: "none", outline: "none",
          transition: "all 0.18s",
          background: waarde === opt.waarde ? "#1a3a5c" : "#f9fafb",
          color: waarde === opt.waarde ? "#fff" : "#4a5568",
          fontFamily: "inherit",
        }}>{opt.label}</button>
      ))}
    </div>
  );
}

function SectieTitel({ icoon, tekst }) {
  return (
    <div style={{ display: "flex", alignItems: "center", gap: 7, margin: "18px 0 11px" }}>
      <span style={{ fontSize: 15 }}>{icoon}</span>
      <span style={{ fontSize: 10, fontWeight: 700, color: "#1a3a5c", letterSpacing: "0.12em", textTransform: "uppercase" }}>
        {tekst}
      </span>
      <div style={{ flex: 1, height: 1, background: "#e2e8f0", marginLeft: 4 }} />
    </div>
  );
}

// ─────────────────────────────────────────────────────────────────────────────
// HOOFDAPPLICATIE
// ─────────────────────────────────────────────────────────────────────────────
const STANDAARD = {
  leeftijd: "50", vrouw: false, dmDuur: "15",
  sbp: "130", ldl: "2.8", hba1c: "68",
  acr: "1.5", egfr: "90",
  roken: false, bewegen: true, eerderCVD: false,
};

// Hulpfunctie: string → getal voor berekeningen
function n(v) { return parseFloat(v); }

export default function App() {
  const [inv, setInv] = useState(STANDAARD);
  const [berekend, setBerekend] = useState(false);
  const [res, setRes] = useState(null);
  const [fouten, setFouten] = useState({});

  const stel = (sleutel) => (e) => {
    // Numerieke velden als string opslaan zodat het veld leeg kan zijn
    const val = e.target
      ? (e.target.type === "number" ? e.target.value : e.target.value)
      : e;
    setInv(p => ({ ...p, [sleutel]: val }));
    setBerekend(false);
  };

  const valideer = () => {
    const f = {};
    const l = n(inv.leeftijd);
    if (isNaN(l) || l < 30 || l > 69)         f.leeftijd = "30–69 jaar vereist";
    const d = n(inv.dmDuur);
    if (isNaN(d) || d < 0 || d > 70)          f.dmDuur = "0–70 jaar";
    const s = n(inv.sbp);
    if (isNaN(s) || s < 80 || s > 220)        f.sbp = "80–220 mmHg";
    const ldl = n(inv.ldl);
    if (isNaN(ldl) || ldl < 0.5 || ldl > 12) f.ldl = "0,5–12 mmol/L";
    const h = n(inv.hba1c);
    if (isNaN(h) || h < 20 || h > 200)        f.hba1c = "20–200 mmol/mol";
    const a = n(inv.acr);
    if (isNaN(a) || a < 0 || a > 500)         f.acr = "0–500 mg/mmol";
    const e = n(inv.egfr);
    if (isNaN(e) || e < 5 || e > 200)         f.egfr = "5–200 mL/min/1,73m²";
    return f;
  };

  // Invoer met numerieke strings omzetten voor berekeningen
  const invNum = () => ({
    ...inv,
    leeftijd: n(inv.leeftijd),
    dmDuur:   n(inv.dmDuur),
    sbp:      n(inv.sbp),
    ldl:      n(inv.ldl),
    hba1c:    n(inv.hba1c),
    acr:      n(inv.acr),
    egfr:     n(inv.egfr),
  });

  const bereken = () => {
    const f = valideer();
    setFouten(f);
    if (Object.keys(f).length > 0) return;
    const i = invNum();
    setRes({
      cvd5:   berekenCVDRisico(i, 5),
      cvd10:  berekenCVDRisico(i, 10),
      cvd30:  berekenCVDRisico(i, 30),
      eskd10: berekenESKDRisico(i, 10),
    });
    setBerekend(true);
  };

  const herstel = () => {
    setInv(STANDAARD);
    setBerekend(false);
    setRes(null);
    setFouten({});
  };

  const acrNum = n(inv.acr);
  const albInfo = albLabel(isNaN(acrNum) ? 0 : acrNum, inv.vrouw);

  return (
    <div style={{
      minHeight: "100vh",
      background: "linear-gradient(145deg, #0b1829 0%, #162d4a 55%, #0b2035 100%)",
      fontFamily: "'DM Sans', 'Helvetica Neue', sans-serif",
      padding: "28px 14px 52px",
    }}>
      <style>{`
        @import url('https://fonts.googleapis.com/css2?family=DM+Serif+Display&family=DM+Sans:opsz,wght@9..40,400;9..40,500;9..40,600;9..40,700&display=swap');
        * { box-sizing: border-box; }
        input[type=number] { -moz-appearance: textfield; }
        input[type=number]::-webkit-inner-spin-button { opacity: 0.4; }
        input:focus, select:focus {
          border-color: #3b82f6 !important;
          background: #fff !important;
          box-shadow: 0 0 0 3px rgba(59,130,246,0.13);
        }
        .kaart {
          background: rgba(255,255,255,0.97);
          border-radius: 18px;
          box-shadow: 0 16px 48px rgba(0,0,0,0.28);
        }
        button { transition: opacity 0.15s, transform 0.1s; }
        button:active { transform: scale(0.98); }
      `}</style>

      {/* Koptekst */}
      <div style={{ textAlign: "center", marginBottom: 24 }}>
        <div style={{ color: "#60a5fa", fontSize: 10, fontWeight: 700, letterSpacing: "0.2em", textTransform: "uppercase", marginBottom: 8 }}>
          Type 1 Diabetes · Steno Risicomodel
        </div>
        <h1 style={{
          fontFamily: "'DM Serif Display', serif",
          color: "#fff", fontSize: "clamp(23px, 5vw, 37px)",
          margin: "0 0 9px", fontWeight: 400, lineHeight: 1.15,
        }}>
          Hart- &amp; Nierziekte<br />Risicocalculator
        </h1>
        <p style={{ color: "#93b8d8", fontSize: 13, margin: "0 auto", maxWidth: 460, lineHeight: 1.6 }}>
          Berekent 5, 10 en 30-jaars HVZ-risico en 10-jaars ESKD-risico
          voor patiënten van 30–69 jaar met type 1 diabetes.
        </p>
      </div>

      <div style={{ maxWidth: 820, margin: "0 auto" }}>

        {/* ── INVOERKAART ── */}
        <div className="kaart" style={{ padding: "22px 22px 18px" }}>

          <SectieTitel icoon="👤" tekst="Patiëntprofiel" />
          <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(148px, 1fr))", gap: 12 }}>
            <Veld label="Leeftijd" hint="(30–69 jr)" fout={fouten.leeftijd}>
              <input type="number" value={inv.leeftijd} min={30} max={69} step={1}
                onChange={stel("leeftijd")}
                style={{ ...invoerStijl, borderColor: fouten.leeftijd ? "#e53935" : "#dde1e7" }} />
            </Veld>
            <Veld label="Geslacht">
              <Schakelaar
                waarde={inv.vrouw ? "vrouw" : "man"}
                onChange={v => { setInv(p => ({ ...p, vrouw: v === "vrouw" })); setBerekend(false); }}
                opties={[{ waarde: "man", label: "Man" }, { waarde: "vrouw", label: "Vrouw" }]} />
            </Veld>
            <Veld label="DM-duur" hint="(jaren)" fout={fouten.dmDuur}>
              <input type="number" value={inv.dmDuur} min={0} max={70} step={1}
                onChange={stel("dmDuur")}
                style={{ ...invoerStijl, borderColor: fouten.dmDuur ? "#e53935" : "#dde1e7" }} />
            </Veld>
          </div>

          <SectieTitel icoon="🔬" tekst="Klinische metingen" />
          <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(148px, 1fr))", gap: 12 }}>
            <Veld label="Systolische BD" hint="(mmHg)" fout={fouten.sbp}>
              <input type="number" value={inv.sbp} step={1}
                onChange={stel("sbp")}
                style={{ ...invoerStijl, borderColor: fouten.sbp ? "#e53935" : "#dde1e7" }} />
            </Veld>
            <Veld label="LDL-cholesterol" hint="(mmol/L)" fout={fouten.ldl}>
              <input type="number" value={inv.ldl} step={0.1}
                onChange={stel("ldl")}
                style={{ ...invoerStijl, borderColor: fouten.ldl ? "#e53935" : "#dde1e7" }} />
            </Veld>
            <Veld label="HbA1c" hint="(mmol/mol)" fout={fouten.hba1c}>
              <input type="number" value={inv.hba1c} step={1}
                onChange={stel("hba1c")}
                style={{ ...invoerStijl, borderColor: fouten.hba1c ? "#e53935" : "#dde1e7" }} />
            </Veld>
            <Veld label="eGFR" hint="(mL/min/1,73m²)" fout={fouten.egfr}>
              <input type="number" value={inv.egfr} step={1}
                onChange={stel("egfr")}
                style={{ ...invoerStijl, borderColor: fouten.egfr ? "#e53935" : "#dde1e7" }} />
            </Veld>
            <Veld label="ACR" hint="(mg/mmol)" fout={fouten.acr}>
              <input type="number" value={inv.acr} step={0.1} min={0}
                onChange={stel("acr")}
                style={{ ...invoerStijl, borderColor: fouten.acr ? "#e53935" : "#dde1e7" }} />
              {!fouten.acr && inv.acr >= 0 && (
                <div style={{ fontSize: 11, marginTop: 4, fontWeight: 600, color: albInfo.kleur }}>
                  → {albInfo.tekst}
                </div>
              )}
            </Veld>
          </div>

          <SectieTitel icoon="⚠️" tekst="Risicofactoren" />
          <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(148px, 1fr))", gap: 12 }}>
            <Veld label="Roken">
              <Schakelaar waarde={inv.roken ? "ja" : "nee"}
                onChange={v => { setInv(p => ({ ...p, roken: v === "ja" })); setBerekend(false); }}
                opties={[{ waarde: "nee", label: "Nee" }, { waarde: "ja", label: "Ja" }]} />
            </Veld>
            <Veld label="Regelmatig bewegen" hint="(≥30 min/dag)">
              <Schakelaar waarde={inv.bewegen ? "ja" : "nee"}
                onChange={v => { setInv(p => ({ ...p, bewegen: v === "ja" })); setBerekend(false); }}
                opties={[{ waarde: "ja", label: "Ja" }, { waarde: "nee", label: "Nee" }]} />
            </Veld>
            <Veld label="Eerder HVZ">
              <Schakelaar waarde={inv.eerderCVD ? "ja" : "nee"}
                onChange={v => { setInv(p => ({ ...p, eerderCVD: v === "ja" })); setBerekend(false); }}
                opties={[{ waarde: "nee", label: "Nee" }, { waarde: "ja", label: "Ja" }]} />
            </Veld>
          </div>

          {/* Knoppen */}
          <div style={{ display: "flex", gap: 10, marginTop: 20 }}>
            <button onClick={bereken} style={{
              flex: 1, padding: "13px 0", borderRadius: 10, border: "none",
              background: "linear-gradient(135deg, #1a3a5c 0%, #2563a8 100%)",
              color: "#fff", fontSize: 15, fontWeight: 700, cursor: "pointer",
              letterSpacing: "0.03em", fontFamily: "inherit",
              boxShadow: "0 4px 18px rgba(37,99,168,0.38)",
            }}>
              Bereken risico
            </button>
            <button onClick={herstel} style={{
              padding: "13px 18px", borderRadius: 10,
              border: "1.5px solid #dde1e7", background: "#f9fafb",
              color: "#4a5568", fontSize: 13, fontWeight: 600,
              cursor: "pointer", fontFamily: "inherit",
            }}>
              Wissen
            </button>
          </div>
        </div>

        {/* ── RESULTATEN ── */}
        {berekend && res && (
          <div style={{ marginTop: 18 }}>

            {/* HVZ */}
            <div className="kaart" style={{ padding: "22px 22px", marginBottom: 14 }}>
              <div style={{ display: "flex", alignItems: "center", gap: 10, marginBottom: 16 }}>
                <div style={{
                  width: 38, height: 38, borderRadius: "50%",
                  background: "linear-gradient(135deg, #b91c1c, #ef4444)",
                  display: "flex", alignItems: "center", justifyContent: "center",
                  fontSize: 19, flexShrink: 0,
                }}>❤️</div>
                <div>
                  <div style={{ fontFamily: "'DM Serif Display', serif", fontSize: 19, color: "#111" }}>
                    Hart- en vaatziekten (HVZ)
                  </div>
                  <div style={{ fontSize: 11, color: "#94a3b8", marginTop: 1 }}>
                    Ischemische hartziekte · Beroerte · Hartfalen · Perifeer arterieel vaatlijden
                  </div>
                </div>
              </div>
              <div style={{ display: "flex", gap: 10, flexWrap: "wrap" }}>
                <Meter pct={res.cvd5}  label="5-jaarsrisico"  sublabel="HVZ" />
                <Meter pct={res.cvd10} label="10-jaarsrisico" sublabel="HVZ" />
                <Meter pct={res.cvd30} label="30-jaarsrisico" sublabel="HVZ" />
              </div>
            </div>

            {/* ESKD */}
            <div className="kaart" style={{ padding: "22px 22px", marginBottom: 14 }}>
              <div style={{ display: "flex", alignItems: "center", gap: 10, marginBottom: 16 }}>
                <div style={{
                  width: 38, height: 38, borderRadius: "50%",
                  background: "linear-gradient(135deg, #1d4ed8, #3b82f6)",
                  display: "flex", alignItems: "center", justifyContent: "center",
                  fontSize: 19, flexShrink: 0,
                }}>🫘</div>
                <div>
                  <div style={{ fontFamily: "'DM Serif Display', serif", fontSize: 19, color: "#111" }}>
                    Eindstadium nierziekte (ESKD)
                  </div>
                  <div style={{ fontSize: 11, color: "#94a3b8", marginTop: 1 }}>
                    Concurrerende risico's model (ESKD versus niet-ESKD sterfte)
                  </div>
                </div>
              </div>
              <div style={{ display: "flex", gap: 10, flexWrap: "wrap", alignItems: "stretch" }}>
                <Meter pct={res.eskd10} label="10-jaarsrisico" sublabel="ESKD" />
                <div style={{
                  flex: 2, minWidth: 200,
                  background: "#f0f5ff", borderRadius: 12, padding: "14px 16px",
                  border: "1px solid #c7d7f8", fontSize: 12, color: "#334155", lineHeight: 1.75,
                }}>
                  <div style={{ fontWeight: 700, color: "#1e3a6e", marginBottom: 5 }}>Methodologie</div>
                  Het ESKD-risico houdt rekening met het concurrerende risico op overlijden zonder ESKD:
                  <div style={{ fontFamily: "monospace", fontSize: 11, color: "#1e3a6e", margin: "6px 0" }}>
                    F<sub>ESKD</sub>(t) = λ<sub>ESKD</sub>/(λ<sub>ESKD</sub>+λ<sub>d</sub>)
                    &nbsp;×&nbsp;[1−e<sup>−(λ<sub>ESKD</sub>+λ<sub>d</sub>)t</sup>]
                  </div>
                  <div style={{ fontSize: 11, color: "#475569", borderTop: "1px solid #dbeafe", paddingTop: 6, marginTop: 4 }}>
                    ACR {inv.acr} mg/mmol →&nbsp;

                    <strong style={{ color: albInfo.kleur }}>{albInfo.tekst}</strong>
                    <span style={{ color: "#94a3b8", marginLeft: 6 }}>
                      (grens: {inv.vrouw ? "3,5 / 35" : "2,5 / 25"} mg/mmol)
                    </span>
                  </div>
                </div>
              </div>
            </div>

            {/* Legenda */}
            <div className="kaart" style={{ padding: "12px 20px" }}>
              <div style={{ display: "flex", gap: 16, flexWrap: "wrap", alignItems: "center" }}>
                <span style={{ fontSize: 10, fontWeight: 700, color: "#94a3b8", textTransform: "uppercase", letterSpacing: "0.1em" }}>
                  Risicocategorieën:
                </span>
                {[
                  { label: "Laag",      bg: "#e8f5e9", text: "#2e7d32", bereik: "<10%" },
                  { label: "Matig",     bg: "#fff8e1", text: "#f57f17", bereik: "10–20%" },
                  { label: "Hoog",      bg: "#fff3e0", text: "#e65100", bereik: "20–40%" },
                  { label: "Zeer hoog", bg: "#fce4ec", text: "#b71c1c", bereik: ">40%" },
                ].map(c => (
                  <div key={c.label} style={{ display: "flex", alignItems: "center", gap: 5 }}>
                    <div style={{ width: 11, height: 11, borderRadius: 3, background: c.bg, border: `1.5px solid ${c.text}` }} />
                    <span style={{ fontSize: 12, color: "#475569" }}><strong>{c.label}</strong> {c.bereik}</span>
                  </div>
                ))}
              </div>
            </div>

            <p style={{ color: "#5a7a96", fontSize: 11, textAlign: "center", marginTop: 14, lineHeight: 1.6, padding: "0 8px" }}>
              Gebaseerd op de Steno Diabetes Center Copenhagen voorspellingsmodellen voor type 1 diabetes (leeftijd 30–69 jaar).
              Uitsluitend bedoeld als klinische beslissingsondersteuning — geen vervanging voor klinisch oordeel.
            </p>
          </div>
        )}
      </div>
    </div>
  );
}
