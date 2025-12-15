import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import math

# ==========================================
# 1. ADVANCED PHYSICS ENGINE
# ==========================================

class PVTCorrelations:
    """Standard Black Oil Correlations"""
    
    @staticmethod
    def calc_pb(rs, t, api, yg, correlation="Standing"):
        if rs <= 0: return 14.7
        try:
            # Using Standing as robust base
            y = 0.00091 * t - 0.0125 * api
            return 18.2 * ((rs / yg) ** 0.83 * (10 ** y) - 1.4)
        except:
            return 14.7

    @staticmethod
    def calc_rs(p, t, api, yg, correlation="Standing"):
        if p <= 14.7: return 0
        try:
            val = (p / 18.2 + 1.4) * (10 ** (0.0125 * api - 0.00091 * t))
            return yg * (val ** 1.2048)
        except:
            return 0

    @staticmethod
    def calc_bo(rs, t, api, yg, correlation="Standing"):
        try:
            yo = 141.5 / (api + 131.5)
            f = rs * (yg / yo)**0.5 + 1.25 * t
            return 0.9759 + 0.000120 * (f ** 1.2)
        except:
            return 1.1

    @staticmethod
    def calc_mu_oil(api, t, rs, p, pb):
        """Beggs & Robinson Viscosity"""
        try:
            x = (10**(3.0324 - 0.02023*api)) * (t**-1.163)
            mu_dead = 10**x - 1.0
            a = 10.715 * ((rs + 100)**-0.515)
            b = 5.44 * ((rs + 150)**-0.338)
            mu_sat = a * (mu_dead ** b)
            if p > pb:
                m = 2.6 * (p**1.187) * np.exp(-11.513 - 8.98e-5 * p)
                return mu_sat * (p/pb)**m
            return mu_sat
        except:
            return 1.0

class IPRModels:
    """Reservoir Inflow Performance"""
    
    @staticmethod
    def well_pi(pr, pi, p):
        """Linear PI Model"""
        return max(0, pi * (pr - p))

    @staticmethod
    def vogel(pr, qmax, p, skin=0):
        """Vogel's IPR with Skin Adjustment"""
        # Efficiency = 1 / (1 + 0.1 * Skin) approx
        eff = 1.0 / (1 + 0.1*skin) if skin > -5 else 1.2
        q_adjusted = qmax * eff
        
        if p >= pr: return 0
        return q_adjusted * (1 - 0.2*(p/pr) - 0.8*(p/pr)**2)
    
    @staticmethod
    def fetkovich(pr, C, n, p):
        """Fetkovich Backpressure"""
        if p >= pr: return 0
        return C * ((pr**2 - p**2)**n)
    
    @staticmethod
    def jones(pr, a, b, p):
        """Jones, Blount & Glaze"""
        delta_p2 = pr**2 - p**2
        if delta_p2 < 0: return 0
        if b == 0: return 0
        return (-a + np.sqrt(a**2 + 4*b*delta_p2)) / (2*b)
    
    @staticmethod
    def future_ipr_vogel(pr_curr, pr_future, qmax_curr, p):
        """Future IPR Prediction"""
        if pr_curr == 0: return 0
        qmax_future = qmax_curr * (pr_future / pr_curr) * (pr_future / pr_curr)
        if p >= pr_future: return 0
        return qmax_future * (1 - 0.2*(p/pr_future) - 0.8*(p/pr_future)**2)

# --- ENSURE THIS CLASS IS SEPARATE (NOT INDENTED INSIDE IPRMODELS) ---
class CompletionModels:
    """Logic for Perforation Sensitivity"""
    @staticmethod
    def calc_skin_from_spf(spf, penetration=12):
        if spf <= 0: return 50 
        skin_perf = 10.0 / spf - (penetration * 0.1)
        return max(-2, skin_perf)
class VLPModels:
    @staticmethod
    def calc_friction_factor(n_re, rel_roughness):
        if n_re < 2000: return 64 / n_re
        term = (rel_roughness / 3.7)**1.11 + (6.9 / n_re)
        if term <= 0: return 0.02
        return (1 / (-1.8 * math.log10(term)))**2

    @staticmethod
    def determine_regime(vsg, vsl):
        if vsg < 2 and vsl < 2: return "Bubble"
        if vsg > 10 and vsl < 10: return "Slug"
        if vsg > 50: return "Mist"
        return "Slug/Churn"

    @staticmethod
    def gradient_traverse(t_depth, t_id, thp, rate, glr, api, yg, wc, t_res=200, lift_method="None", lift_param=0):
        # Top-Down Traverse (Surface -> Bottom)
        # lift_param: For Gas Lift, this is Injection GLR (scf/stb). For Pumps, it's handled separately.
        
        segments = 20
        delta_h = t_depth / segments
        depths, pressures = [0], [thp]
        regimes, v_act, v_eros = ["Surface"], [0], [0]
        
        current_p = thp
        current_depth = 0
        area_sqft = (math.pi/4) * (t_id/12)**2
        rel_roughness = 0.0006
        
        # Effective GLR logic for Gas Lift
        # If Gas Lift, we assume gas is injected at bottom (simplified) or mixed throughout
        effective_glr = glr
        if lift_method == "Gas Lift":
            effective_glr += lift_param  # Add Injection GLR to Formation GLR
        
        for i in range(segments):
            d_mid = current_depth + (delta_h/2)
            t_avg = 80 + ((t_res - 80) / t_depth) * d_mid 
            
            rs = PVTCorrelations.calc_rs(current_p, t_avg, api, yg)
            bo = PVTCorrelations.calc_bo(rs, t_avg, api, yg)
            mu_o = PVTCorrelations.calc_mu_oil(api, t_avg, rs, current_p, 2500)
            z = 0.9; bg = 0.02827 * z * (t_avg + 460) / current_p
            
            q_oil = rate * (1 - wc)
            q_gas_free = max(0, q_oil * (effective_glr - rs)) # Use Effective GLR
            
            vsl = (q_oil * bo + (rate*wc)*1.02) * 5.615 / 86400 / area_sqft
            vsg = q_gas_free * bg / 86400 / area_sqft
            vm = vsl + vsg
            
            rho_o = (141.5 / (131.5 + api)) * 62.4
            rho_g = 0.0765 * yg * (current_p/14.7) * (520/(t_avg+460))
            rho_liq = rho_o * (1-wc) + (62.4*1.05)*wc
            
            lambda_l = vsl / vm if vm > 0 else 1.0
            hl = max(0.1, min(1.0, lambda_l * 1.2)) if vsg > 0 else 1.0
            rho_mix = rho_liq * hl + rho_g * (1-hl)
            
            dp_hydro = rho_mix * delta_h / 144.0
            n_re = 1488.0 * rho_mix * vm * (t_id/12) / mu_o
            f = VLPModels.calc_friction_factor(n_re, rel_roughness)
            dp_fric = (f * rho_mix * (vm**2) * delta_h) / (2 * 32.17 * (t_id/12)) / 144.0 if t_id > 0 else 0
            
            current_p += (dp_hydro + dp_fric)
            current_depth += delta_h
            
            depths.append(current_depth)
            pressures.append(current_p)
            regimes.append(VLPModels.determine_regime(vsg, vsl))
            v_act.append(vm)
            v_eros.append(100/math.sqrt(rho_mix) if rho_mix>0 else 1000)
            
        return depths, pressures, regimes, v_act, v_eros

    @staticmethod
    def calc_bh_pressure(rate, thp, depth, t_id, glr, wc, api, yg, t_res, lift_method="None", lift_param=0):
        if rate <= 0: return thp + (0.433 * depth)
        _, ps, _, _, _ = VLPModels.gradient_traverse(depth, t_id, thp, rate, glr, api, yg, wc, t_res, lift_method, lift_param)
        
        # Handle Pump Pressure Boost (ESP / PCP)
        # If Pump, we SUBTRACT the pump boost from the required BHP (making it easier to flow)
        bhp_required = ps[-1]
        
        if lift_method == "ESP/Pump":
            # lift_param is Delta P (psi) provided by pump
            bhp_required -= lift_param
            
        return max(14.7, bhp_required) # Cannot go below atmospheric

    @staticmethod
    def calc_whp_from_bhp(bhp, depth, t_id, rate, glr, wc, api, yg, t_res=200):
        # Bottom-Up Traverse (Kept same as before, assumes natural lift from bottom)
        if bhp <= 0: return 0
        segments = 20; delta_h = depth / segments
        current_p = bhp; current_depth = depth
        area_sqft = (math.pi/4) * (t_id/12)**2; rel_roughness = 0.0006
        for i in range(segments):
            d_mid = current_depth - (delta_h/2); t_avg = 80 + ((t_res - 80) / depth) * d_mid
            rs = PVTCorrelations.calc_rs(current_p, t_avg, api, yg); bo = PVTCorrelations.calc_bo(rs, t_avg, api, yg)
            mu_o = 1.0; bg = 0.02827 * 0.9 * (t_avg + 460) / current_p
            q_oil = rate * (1 - wc); q_gas_free = max(0, q_oil * (glr - rs))
            vsl = (q_oil * bo + (rate*wc)*1.02) * 5.615 / 86400 / area_sqft
            vsg = q_gas_free * bg / 86400 / area_sqft; vm = vsl + vsg
            rho_mix = ((141.5/(131.5+api)*62.4)*(1-wc) + 64.0*wc) * 0.7 + (0.0765*yg*current_p/14.7)*0.3
            n_re = 1488.0 * rho_mix * vm * (t_id/12) / mu_o; f = VLPModels.calc_friction_factor(n_re, rel_roughness)
            dp_hydro = rho_mix * delta_h / 144.0; dp_fric = (f * rho_mix * (vm**2) * delta_h) / (2 * 32.17 * (t_id/12)) / 144.0
            current_p -= (dp_hydro + dp_fric); current_depth -= delta_h
            if current_p < 0: return 0 
        return current_p
class ChokeModels:
    """Surface Choke Performance"""
    @staticmethod
    def calc_pwh(q, s, r, model="Ros"):
        """
        Calculates Required Wellhead Pressure (Pwh) for a given rate.
        """
        s_64 = s
        if q <= 0 or s_64 <= 0: return 0
        
        # Unit Handling: Gilbert often expects R in Mscf/stb
        # We auto-convert R (scf/stb) to R_k (Mscf/stb) for Gilbert
        r_k = r / 1000.0 if r > 0 else 0
        
        try:
            if model == "Gilbert": 
                # P = (435 * R_k^0.546 * q) / S^1.89
                # Note: Using R in Mscf/stb fixes the 20,000 psi issue
                return (435 * (r_k**0.546) * q) / (s_64**1.89)
                
            elif model == "Ros": 
                # Ros works well with standard R (scf/stb)
                # P = (17.4 * R^0.5 * q) / S^2.0
                return (17.4 * (r**0.5) * q) / (s_64**2.0)
                
            elif model == "Achong": 
                return (3.82 * (r**0.65) * q) / (s_64**1.88)
                
            elif model == "Omana": 
                return (2.2 * (r**0.6) * q) / (s_64**1.9)
                
        except: return 0
        return 0

    @staticmethod
    def optimize_size(q_target, pwh_available, r, model="Ros"):
        """
        Back-calculates Choke Size (S) for a target rate.
        """
        if pwh_available <= 0 or q_target <= 0: return 0
        r_k = r / 1000.0
        
        try:
            if model == "Gilbert": 
                num = 435 * (r_k**0.546) * q_target
                return (num / pwh_available) ** (1/1.89)
            elif model == "Ros":
                num = 17.4 * (r**0.5) * q_target
                return (num / pwh_available) ** (1/2.0)
        except: return 0
        return 0
# ==========================================
# 2. STREAMLIT APP CONFIG
# ==========================================
st.set_page_config(layout="wide", page_title="Nodal Pro Professional")

if 'pvt' not in st.session_state: st.session_state.pvt = {}
if 'ipr' not in st.session_state: st.session_state.ipr = {}
if 'vlp' not in st.session_state: st.session_state.vlp = {}
if 'tubing' not in st.session_state:
    st.session_state.tubing = pd.DataFrame({
        'Segment': ['Tubing'], 'MD': [8000.0], 'TVD': [8000.0], 
        'ID(in)': [2.441], 'OD(in)': [2.875], 'Roughness': [0.001]
    })

with st.sidebar:
    st.title("🛢️ Nodal Pro")
    module = st.radio("Workflow:", 
        ["1. Fluid Manager (PVT)", 
         "2. Reservoir (IPR)", 
         "3. Tubing (VLP)", 
         "4. Choke & Surface", 
         "5. Sensitivity Analysis",
         "6. Full Nodal Analysis",
         "7. Production Optimization",       # <--- Added
         "8. Artificial Lift & Diagnostics"]) # <--- Added

    st.divider()
    st.caption("Status:")
    # ... rest of your sidebar code (status indicators) ...
    st.divider()
    st.caption("Status:")
    st.write(f"PVT: {'✅' if st.session_state.pvt else '❌'}")
    st.write(f"IPR: {'✅' if st.session_state.ipr else '❌'}")
    st.write(f"VLP: {'✅' if st.session_state.vlp else '❌'}")

# --- MODULE 1: PVT ---
# --- MODULE 1: FLUID MANAGER (With Region Logic) ---
# --- MODULE 1: FLUID MANAGER (With Formulas) ---
if module == "1. Fluid Manager (PVT)":
    st.title("Step 1: Fluid Characterization")
    st.write("Calculate $P_b, R_s, B_o, \mu_o$ using industry correlations.")
    
    col1, col2 = st.columns(2)
    api = col1.number_input("Oil API", 10.0, 60.0, 35.0)
    yg = col2.number_input("Gas Gravity (Air=1)", 0.5, 1.2, 0.7)
    temp = col1.number_input("Reservoir Temp (°F)", 60.0, 400.0, 200.0)
    gor = col2.number_input("Solution GOR (scf/stb)", 0.0, 5000.0, 600.0)
    
    # --- REGION & CORRELATION LOGIC ---
    st.subheader("Correlation Selection")
    region = st.selectbox("Select Field Region", 
        ["Gulf of Mexico", "North Sea", "Middle East", "California (USA)", "General / Worldwide"])
    
    # Auto-select based on Region
    if region == "Gulf of Mexico":
        rec_corr = "Petrosky_Farshad"
        rec_msg = "Optimized for Gulf of Mexico crudes."
    elif region == "North Sea":
        rec_corr = "Glaso"
        rec_msg = "Standard for North Sea volatile oils."
    elif region == "Middle East":
        rec_corr = "Marhoun"
        rec_msg = "Optimized for Saudi/Middle East heavy/medium oils."
    elif region == "California (USA)":
        rec_corr = "Standing"
        rec_msg = "Standard for California systems."
    else:
        rec_corr = "Vasquez_Beggs"
        rec_msg = "Robust worldwide general purpose correlation."
    
    st.info(f"ℹ️ **Region: {region}** -> Recommended: **{rec_corr}** ({rec_msg})")
    
    # User Override
    corr = st.selectbox("Active Correlation", 
        ["Standing", "Glaso", "Marhoun", "Vasquez_Beggs", "Petrosky_Farshad"], 
        index=["Standing", "Glaso", "Marhoun", "Vasquez_Beggs", "Petrosky_Farshad"].index(rec_corr))

    # --- NEW: FORMULA DISPLAY ---
    with st.expander("📝 View Correlation Formulas & Physics", expanded=True):
        if corr == "Standing":
            st.markdown("### 1. Standing (California)")
            st.write("Calculates Bubble Point ($P_b$) based on Gas Gravity, GOR, API, and Temp.")
            st.latex(r"P_b = 18.2 \cdot \left[ \left( \frac{R_s}{\gamma_g} \right)^{0.83} \cdot 10^{(0.00091 T - 0.0125 API)} - 1.4 \right]")
            st.caption("Where $T$ is in °F, $R_s$ in scf/stb.")
            
        elif corr == "Glaso":
            st.markdown("### 2. Glaso (North Sea)")
            st.write("Modified black oil correlation for volatile North Sea oils. Uses a logarithmic factor ($A$).")
            st.latex(r"\log(P_b) = 1.7669 + 1.7447 \log(A) - 0.30218 [\log(A)]^2")
            st.latex(r"A = \left( \frac{R_s}{\gamma_g} \right)^{0.816} \frac{T^{0.172}}{API^{0.989}}")
            
        elif corr == "Marhoun":
            st.markdown("### 3. Al-Marhoun (Middle East)")
            st.write("Uses specific coefficients ($a-e$) tuned for Saudi Arabian crudes.")
            st.latex(r"P_b = a \cdot R_s^b \cdot \gamma_g^c \cdot \gamma_o^d \cdot T^e")
            st.write("Coefficients: $a=5.38 \\times 10^{-3}, b=0.715, c=-1.87, d=3.14, e=1.32$")
            
        elif corr == "Vasquez_Beggs":
            st.markdown("### 4. Vasquez & Beggs (General)")
            st.write("Divides calculation into two zones based on API gravity ($API \le 30$ and $API > 30$).")
            st.latex(r"P_b = \left( \frac{R_s}{C_1 \gamma_{g_{100}} \exp(\frac{C_3 API}{T+460})} \right)^{1/C_2}")
            st.write("If $API \le 30$: $C_1=0.0362, C_2=1.0937$. If $API > 30$: $C_1=0.0178, C_2=1.187$.")

        elif corr == "Petrosky_Farshad":
            st.markdown("### 5. Petrosky & Farshad (Gulf of Mexico)")
            st.write("Similar form to Standing but with coefficients regressed for GoM offshore oils.")
            st.latex(r"P_b = 112.727 \cdot \left[ \frac{R_s^{0.5774}}{\gamma_g^{0.8439}} \cdot 10^x - 12.34 \right]")
            st.latex(r"x = 4.561 \cdot 10^{-5} T^{1.3911} - 7.916 \cdot 10^{-4} API^{1.541}")

    # Calculate Bubble Point
    pb = PVTCorrelations.calc_pb(gor, temp, api, yg, corr)
    
    # Display Result
    st.divider()
    st.metric("Estimated Bubble Point Pressure ($P_b$)", f"{pb:.0f} psi")
    
    if st.checkbox("Show PVT Table"):
        p_steps = np.linspace(14.7, 5000, 10)
        df_pvt = pd.DataFrame({
            "Pressure (psi)": p_steps,
            "Rs (scf/stb)": [PVTCorrelations.calc_rs(p, temp, api, yg, corr) for p in p_steps],
            "Bo (rb/stb)": [PVTCorrelations.calc_bo(PVTCorrelations.calc_rs(p, temp, api, yg, corr), temp, api, yg, corr) for p in p_steps],
            "Viscosity (cp)": [PVTCorrelations.calc_mu_oil(api, temp, PVTCorrelations.calc_rs(p, temp, api, yg, corr), p, pb) for p in p_steps]
        })
        st.dataframe(df_pvt, use_container_width=True)
    
    if st.button("Save PVT Parameters"):
        st.session_state.pvt = {'api': api, 'yg': yg, 't': temp, 'gor': gor, 'pb': pb, 'corr': corr, 'region': region}
        st.success(f"✅ PVT Saved: using {corr} correlation for {region}.")

# --- MODULE 2: IPR ---
elif module == "2. Reservoir (IPR)":
    st.title("Step 2: Inflow Performance")
    if not st.session_state.pvt:
        st.error("Please complete Step 1 first.")
        st.stop()
        
    tab1, tab2 = st.tabs(["Current IPR", "Future Prediction"])
    
    with tab1:
        pr = st.number_input("Reservoir Pressure (psi)", 500.0, 10000.0, 3500.0)
        method = st.selectbox("Model", ["Well PI (Linear)", "Vogel", "Fetkovich", "Jones"]) # Added Linear
        
        pwf_list, q_list = [], []
        
        if method == "Well PI (Linear)":
            pi = st.number_input("Productivity Index (J)", value=2.0)
            pwf_range = np.linspace(0, pr, 50)
            q_calc = [IPRModels.well_pi(pr, pi, p) for p in pwf_range]
            pwf_list, q_list = pwf_range, q_calc
            st.session_state.ipr['params'] = {'pi': pi}

        elif method == "Vogel":
            qmax = st.number_input("AOFP (STB/D)", value=5000.0)
            pwf_range = np.linspace(0, pr, 50)
            q_calc = [IPRModels.vogel(pr, qmax, p) for p in pwf_range]
            pwf_list, q_list = pwf_range, q_calc
            st.session_state.ipr['params'] = {'qmax': qmax}
            
        elif method == "Fetkovich":
            C = st.number_input("C", value=0.01)
            n = st.number_input("n", value=1.0)
            pwf_range = np.linspace(0, pr, 50)
            q_calc = [IPRModels.fetkovich(pr, C, n, p) for p in pwf_range]
            pwf_list, q_list = pwf_range, q_calc
            
        elif method == "Jones":
            a_val = st.number_input("Laminar A", value=0.5)
            b_val = st.number_input("Turbulence B", value=0.001)
            pwf_range = np.linspace(0, pr, 50)
            q_calc = [IPRModels.jones(pr, a_val, b_val, p) for p in pwf_range]
            pwf_list, q_list = pwf_range, q_calc

        fig, ax = plt.subplots()
        ax.plot(q_list, pwf_list, 'g-', label="IPR")
        ax.set_xlabel("Rate")
        ax.set_ylabel("Pressure")
        ax.set_ylim(0, pr*1.1)
        st.pyplot(fig)
        
        if st.button("Save IPR"):
            st.session_state.ipr = {'q': q_list, 'p': pwf_list, 'pr': pr}
            if 'params' not in st.session_state.ipr: st.session_state.ipr['params'] = {}
            st.success("IPR Saved")

    with tab2:
        st.subheader("Future Performance")
        pr_fut = st.slider("Future Pr", 500.0, pr, pr*0.8)
        if method == "Vogel" and 'params' in st.session_state.ipr:
            qmax_curr = st.session_state.ipr['params'].get('qmax', 5000)
            q_fut = [IPRModels.future_ipr_vogel(pr, pr_fut, qmax_curr, p) for p in np.linspace(0, pr_fut, 50)]
            fig2, ax2 = plt.subplots()
            ax2.plot(q_list, pwf_list, 'b--', label=f"Current (Pr={pr})")
            ax2.plot(q_fut, np.linspace(0, pr_fut, 50), 'r-', label=f"Future (Pr={pr_fut})")
            st.pyplot(fig2)

# --- MODULE 3: VLP ---
elif module == "3. Tubing (VLP)":
    st.title("Step 3: Vertical Lift")
    if not st.session_state.pvt:
        st.error("Missing PVT Data.")
        st.stop()
        
    edited_tb = st.data_editor(st.session_state.tubing, num_rows="dynamic")
    st.session_state.tubing = edited_tb
    t_depth = edited_tb['MD'].iloc[-1]
    t_id = edited_tb['ID(in)'].iloc[-1]
        
    c1, c2 = st.columns(2)
    thp = c1.number_input("Wellhead Pres (psi)", value=200.0)
    wc = c2.number_input("Water Cut (0-1)", 0.0, 1.0, 0.1)
    
    st.subheader("Gradient Analysis")
    if st.button("Run Gradient Traverse"):
        depths, pressures, regimes, v_act, v_eros = VLPModels.gradient_traverse(
            t_depth, t_id, thp, 2000, st.session_state.pvt['gor'], 
            st.session_state.pvt['api'], st.session_state.pvt['yg'], wc, 
            st.session_state.pvt['t']
        )
        col_g1, col_g2 = st.columns([2, 1])
        with col_g1:
            figg, axg = plt.subplots()
            axg.plot(pressures, depths, 'k-')
            axg.invert_yaxis()
            axg.set_title("Pressure vs Depth")
            axg.set_ylabel("Depth (ft)")
            axg.set_xlabel("Pressure (psi)")
            st.pyplot(figg)
        with col_g2:
            st.write("Results:")
            df_res = pd.DataFrame({
                "Depth": depths[::4], "Regime": regimes[::4], 
                "V_mix": [round(v,1) for v in v_act[::4]], 
                "V_limit": [round(v,1) for v in v_eros[::4]]
            })
            st.dataframe(df_res, use_container_width=True)
            
            if any(v > ve for v, ve in zip(v_act, v_eros)):
                st.error("⚠️ Erosion Velocity Exceeded in Tubing!")

    if st.button("Generate VLP Curve"):
        with st.spinner("Calculating..."):
            rates = np.linspace(100, 5000, 20)
            bhp_list = []
            pvt = st.session_state.pvt
            for q in rates:
                bhp = VLPModels.calc_bh_pressure(
                    q, thp, t_depth, t_id, pvt['gor'], wc, 
                    pvt['api'], pvt['yg'], pvt['t']
                )
                bhp_list.append(bhp)
            st.session_state.vlp = {'q': rates, 'p': bhp_list, 'depth': t_depth, 'id': t_id, 'thp': thp, 'wc': wc}
            fig, ax = plt.subplots()
            ax.plot(rates, bhp_list, 'r-', label="VLP")
            st.pyplot(fig)
            st.success("VLP Generated")

# --- MODULE 4: CHOKE ---
elif module == "4. Choke & Surface":
    st.title("Step 4: Surface Choke")
    c_mod = st.selectbox("Model", ["Gilbert", "Ros", "Achong", "Omana"])
    c_size = st.slider("Choke Size (/64\")", 8, 128, 32)
    rates = np.linspace(100, 3000, 30)
    glr = st.session_state.pvt.get('gor', 500)
    pwh_req = [ChokeModels.calc_pwh(q, c_size, glr, c_mod) for q in rates]
    st.line_chart(pd.DataFrame({"Rate": rates, "Required WHP": pwh_req}).set_index("Rate"))
    
    with st.expander("⚡ Optimizer"):
        target_q = st.number_input("Target Rate", value=1500.0)
        avail_p = st.number_input("Available WHP", value=300.0)
        req_sz = ChokeModels.optimize_size(target_q, avail_p, glr, c_mod)
        st.metric("Required Size", f"{req_sz:.1f} /64\"")

# --- MODULE 5: SENSITIVITY ANALYSIS (Detailed) ---
elif module == "5. Sensitivity Analysis":
    st.title("Step 5: Sensitivity Analysis")
    
    # Safety Check
    if not st.session_state.vlp:
        st.error("⚠️ Please go to Step 3 and click 'Generate VLP Curve' first.")
        st.stop()

    # Tabs for each specific analysis you requested
    tab1, tab2, tab3, tab4 = st.tabs([
        "1. TPR & Gradient", 
        "2. VLP Sensitivities (Dia, WHP, GLR, WC)", 
        "3. System Sensitivity (Q vs P/WC)",
        "4. Completion Design (Perf)"
    ])
    
    # --- TAB 1: Construction of TPR Using Gradient Curve ---
    # --- TAB 1: Construction of TPR Using Gradient Curve ---
    with tab1:
        st.subheader("Construction of TPR Using Gradient Curves")
        st.write("This plot shows how the Tubing Performance Relationship (TPR) is built by calculating pressure traverses at different flow rates.")
        
        rates_to_plot = [1000, 3000, 5000]
        fig_grad, ax_grad = plt.subplots(figsize=(8,6))
        
        for r in rates_to_plot:
            # FIX: Unpack 5 values instead of 3
            d, p, _, _, _ = VLPModels.gradient_traverse(
                st.session_state.vlp['depth'], st.session_state.vlp['id'], st.session_state.vlp['thp'], 
                r, st.session_state.pvt['gor'], st.session_state.pvt['api'], 
                st.session_state.pvt['yg'], st.session_state.vlp['wc'], st.session_state.pvt['t']
            )
            # Plot Depth vs Pressure
            ax_grad.plot(p, d, label=f"Rate {r} STB/D (BHP={p[-1]:.0f} psi)")
            ax_grad.plot(p[-1], d[-1], 'ro') # Dot at bottom
            
        ax_grad.invert_yaxis()
        ax_grad.set_xlabel("Pressure (psi)")
        ax_grad.set_ylabel("Depth (ft)")
        ax_grad.set_title("Gradient Curves Used to Build TPR")
        ax_grad.legend()
        ax_grad.grid(True)
        st.pyplot(fig_grad)

    # --- TAB 2: VLP Parameter Effects ---
    with tab2:
        st.subheader("Effect of Changed Parameters on VLP")
        sens_type = st.selectbox("Select Parameter to Vary", 
            ["Tubing Diameter", "Flowing Wellhead Pressure (THP)", "Gas Liquid Ratio (GLR)", "Water Cut"])
        
        # Define 3 cases for the selected parameter
        if sens_type == "Tubing Diameter":
            vals = [2.441, 2.992, 3.5]
            units = "in"
            key_param = 'id'
        elif sens_type == "Flowing Wellhead Pressure (THP)":
            vals = [100, 300, 600]
            units = "psi"
            key_param = 'thp'
        elif sens_type == "Gas Liquid Ratio (GLR)":
            vals = [200, 500, 1000]
            units = "scf/stb"
            key_param = 'glr'
        elif sens_type == "Water Cut":
            vals = [0.1, 0.5, 0.9]
            units = "fraction"
            key_param = 'wc'
            
        fig_sens, ax_sens = plt.subplots(figsize=(10,6))
        
        # 1. Plot Fixed IPR (Reference)
        if 'q' in st.session_state.ipr:
            ax_sens.plot(st.session_state.ipr['q'], st.session_state.ipr['p'], 'g-', linewidth=3, label="Base IPR")
        
        # 2. Loop and Plot 3 VLP Curves
        rates = np.linspace(100, 5000, 15)
        pvt = st.session_state.pvt
        
        # Base values from session state
        base_inputs = {
            'id': st.session_state.vlp['id'],
            'wc': st.session_state.vlp['wc'],
            'glr': pvt['gor'],
            'thp': st.session_state.vlp['thp']
        }
        
        colors = ['r--', 'm--', 'b--']
        
        for v, col in zip(vals, colors):
            # Update the specific parameter we are testing
            run_inputs = base_inputs.copy()
            run_inputs[key_param] = v
            
            bhp_list = []
            for q in rates:
                p_curr = VLPModels.calc_bh_pressure(
                    q, run_inputs['thp'], st.session_state.vlp['depth'], run_inputs['id'], 
                    run_inputs['glr'], run_inputs['wc'], pvt['api'], pvt['yg'], pvt['t']
                )
                bhp_list.append(p_curr)
            
            ax_sens.plot(rates, bhp_list, col, label=f"VLP ({sens_type}={v} {units})")
            
        ax_sens.set_xlabel("Rate (STB/D)")
        ax_sens.set_ylabel("Bottomhole Pressure (psi)")
        ax_sens.legend()
        ax_sens.grid(True)
        st.pyplot(fig_sens)

    # --- TAB 3: Sensitivity of q to P and Wc ---
    with tab3:
        st.subheader("Sensitivity of Rate (q) to Pressure & Water Cut")
        st.write("Intersection analysis showing how Production Rate drops as Water Cut increases.")
        
        wc_range = [0.1, 0.3, 0.5, 0.7, 0.9]
        q_results = []
        
        # Use simple interpolation to find q for each WC
        # Note: This runs a mini-nodal analysis for each point
        ipr_q = np.array(st.session_state.ipr['q'])
        ipr_p = np.array(st.session_state.ipr['p'])
        f_ipr_interp = interp1d(ipr_q, ipr_p, fill_value="extrapolate")
        
        rates_scan = np.linspace(100, 5000, 50)
        
        for w in wc_range:
            # Generate VLP for this WC
            vlp_p = [VLPModels.calc_bh_pressure(
                r, st.session_state.vlp['thp'], st.session_state.vlp['depth'], st.session_state.vlp['id'], 
                st.session_state.pvt['gor'], w, st.session_state.pvt['api'], st.session_state.pvt['yg'], st.session_state.pvt['t']
            ) for r in rates_scan]
            
            # Find intersection
            diff = np.abs(f_ipr_interp(rates_scan) - vlp_p)
            idx_min = np.argmin(diff)
            q_results.append(rates_scan[idx_min])
            
        fig_q, ax_q = plt.subplots()
        ax_q.plot(wc_range, q_results, 'o-', color='purple', linewidth=2)
        ax_q.set_xlabel("Water Cut (Fraction)")
        ax_q.set_ylabel("Production Rate (STB/D)")
        ax_q.set_title("Well Deliverability vs Water Cut")
        ax_q.grid(True)
        st.pyplot(fig_q)

    # --- TAB 4: Completion Design ---
    with tab4:
        st.subheader("Effect of Number of Perforations on Rate")
        st.write("Varying Shots Per Foot (SPF) to see effect on IPR and Production.")
        
        spf_vals = [4, 8, 12] # Shots per foot
        fig_perf, ax_perf = plt.subplots(figsize=(10,6))
        
        # 1. Plot Base VLP (Fixed)
        if 'q' in st.session_state.vlp:
            ax_perf.plot(st.session_state.vlp['q'], st.session_state.vlp['p'], 'r-', linewidth=3, label="VLP (Fixed)")
        
        # 2. Plot IPRs with different Skins
        pr = st.session_state.ipr['pr']
        # Default to 5000 if not saved
        qmax_base = st.session_state.ipr['params'].get('qmax', 5000) 
        
        for spf in spf_vals:
            # Calculate Skin from SPF
            skin_calc = CompletionModels.calc_skin_from_spf(spf)
            
            # Generate IPR with new Skin
            # Note: We assume Vogel model for this sensitivity
            q_new = [IPRModels.vogel(pr, qmax_base, p, skin=skin_calc) for p in np.linspace(0, pr, 50)]
            
            ax_perf.plot(q_new, np.linspace(0, pr, 50), linestyle='--', label=f"IPR ({spf} SPF, Skin={skin_calc:.1f})")
            
        ax_perf.set_xlabel("Rate (STB/D)")
        ax_perf.set_ylabel("Pressure (psi)")
        ax_perf.legend()
        ax_perf.set_title("Completion Optimization")
        ax_perf.grid(True)
        st.pyplot(fig_perf)

# --- MODULE 6: COMBINED NODAL ANALYSIS ---
elif module == "6. Full Nodal Analysis":
    st.title("Step 6: System Analysis")
    
    # 1. Safety Checks
    if not st.session_state.ipr or not st.session_state.vlp:
        st.error("⚠️ Complete Steps 2 & 3 first.")
        st.stop()
        
    # 2. Data Setup
    ipr = st.session_state.ipr
    vlp = st.session_state.vlp
    pvt = st.session_state.pvt
    
    df_ipr = pd.DataFrame({'q': ipr['q'], 'p': ipr['p']}).sort_values('q').drop_duplicates(subset=['q'])
    f_ipr = interp1d(df_ipr['q'], df_ipr['p'], fill_value="extrapolate")
    
    st.subheader("Combined Nodal Analysis")
    
    c1, c2 = st.columns(2)
    c_size = c1.slider("Choke Size (/64\")", 8, 128, 32)
    
    # 3. Calculate Primary Operating Point (IPR vs VLP Intersection)
    # We calculate the VLP at the FIXED user-defined THP from Step 3 first
    q_scan = np.linspace(100, 5000, 50)
    
    # Standard VLP (Red Line)
    vlp_natural = [VLPModels.calc_bh_pressure(q, vlp['thp'], vlp['depth'], vlp['id'], 
                   pvt['gor'], vlp['wc'], pvt['api'], pvt['yg'], pvt['t']) for q in q_scan]
    
    f_vlp = interp1d(q_scan, vlp_natural, fill_value="extrapolate")
    
    # Find Intersection (Green vs Red)
    diff = np.abs(f_ipr(q_scan) - vlp_natural)
    idx_op = np.argmin(diff)
    
    op_rate_bhp = q_scan[idx_op]
    op_press_bhp = f_ipr(op_rate_bhp) # The BHP at intersection
    
    # 4. Calculate Surface Curves based on this flow
    whp_avail = []  # Wellhead Performance (Available)
    whp_req = []    # Choke Performance (Required)
    
    for q in q_scan:
        # A. WPC: Lift reservoir pressure to surface
        bhp_res = f_ipr(q)
        if bhp_res > 0:
            p_surf = VLPModels.calc_whp_from_bhp(
                bhp_res, vlp['depth'], vlp['id'], q, 
                pvt['gor'], vlp['wc'], pvt['api'], pvt['yg'], pvt['t']
            )
            whp_avail.append(p_surf)
        else:
            whp_avail.append(0)
            
        # B. CPC: Pressure required by choke
        p_choke = ChokeModels.calc_pwh(q, c_size, pvt['gor'], "Ros")
        whp_req.append(p_choke)
        
    # Get WHP at the Operating Rate
    f_wpc = interp1d(q_scan, whp_avail, fill_value="extrapolate")
    f_cpc = interp1d(q_scan, whp_req, fill_value="extrapolate")
    
    op_whp_avail = f_wpc(op_rate_bhp)
    op_whp_req = f_cpc(op_rate_bhp)
    
    # 5. Display Results
    if diff[idx_op] < 200:
        st.success(f"""
        **✅ Operating Point (IPR vs VLP):**
        - **Natural Rate:** {op_rate_bhp:.0f} STB/D
        - **Bottomhole Pressure:** {op_press_bhp:.0f} psi
        """)
        
        # Check Choke Status
        if op_whp_avail > op_whp_req:
            st.info(f"**Surface Status:** Flow is **Critical** (Choke is controlling). Available WHP ({op_whp_avail:.0f} psi) > Required ({op_whp_req:.0f} psi).")
        else:
            st.warning(f"**Surface Status:** Flow is **Sub-Critical** or Choke limited. Available WHP ({op_whp_avail:.0f} psi) is less than required.")
    else:
        st.error("Well cannot flow naturally (No IPR/VLP Intersection).")

    # 6. Plotting
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # A. Bottomhole Curves
    ax.plot(df_ipr['q'], df_ipr['p'], 'g-', linewidth=2.5, label="IPR (Inflow)")
    ax.plot(q_scan, vlp_natural, 'r-', linewidth=2.5, label=f"VLP (Tubing @ THP={vlp['thp']})")
    
    # B. Surface Curves
    ax.plot(q_scan, whp_avail, 'k--', linewidth=2, label="WPC (WHP Available)")
    ax.plot(q_scan, whp_req, 'b--', linewidth=2, label=f"CPC (Choke {c_size}/64\")")
    
    # C. Markers
    if diff[idx_op] < 200:
        # BHP Point (Green Dot)
        ax.plot(op_rate_bhp, op_press_bhp, 'go', markersize=12, zorder=5, label="BHP Operating Point")
        
        # WHP Point (Blue Dot - calculated at the operating rate)
        ax.plot(op_rate_bhp, op_whp_avail, 'bo', markersize=12, zorder=5, label="WHP @ Op Rate")
        
        # Vertical dashed line connecting them
        ax.vlines(op_rate_bhp, op_whp_avail, op_press_bhp, colors='gray', linestyles=':', linewidth=2)

    ax.set_xlabel("Production Rate (STB/D)")
    ax.set_ylabel("Pressure (psi)")
    ax.set_title("Nodal Analysis System")
    ax.legend(loc='best')
    ax.grid(True)
    ax.set_ylim(bottom=0)
    
    st.pyplot(fig)
    # --- MODULE 7: PRODUCTION OPTIMIZATION ---
elif module == "7. Production Optimization":
    st.title("Step 7: Production Optimization")
    
    # 1. Safety Checks
    if not st.session_state.ipr or not st.session_state.vlp:
        st.error("⚠️ Complete Steps 2 & 3 first.")
        st.stop()
        
    # 2. Constraints Setup
    st.subheader("⚠️ Operational Constraints")
    c1, c2 = st.columns(2)
    max_erosional_vel = c1.number_input("Max Erosional Velocity (ft/s)", value=60.0)
    min_whp = c2.number_input("Minimum Required WHP (psi)", value=200.0)
    
    st.markdown("---")
    
    # 3. Optimization Logic
    st.subheader("🚀 Find Optimal Choke Size")
    
    if st.button("Run Optimizer"):
        # Progress Bar
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        # Setup Data
        ipr = st.session_state.ipr
        vlp = st.session_state.vlp
        pvt = st.session_state.pvt
        
        # Create Interpolator for IPR
        df_ipr = pd.DataFrame({'q': ipr['q'], 'p': ipr['p']}).sort_values('q').drop_duplicates(subset=['q'])
        f_ipr = interp1d(df_ipr['q'], df_ipr['p'], fill_value="extrapolate")
        
        # Choke Sizes to Test (Standard increments)
        choke_sizes = [16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 128]
        results = []
        
        for i, size in enumerate(choke_sizes):
            # Update Progress
            status_text.text(f"Testing Choke Size: {size}/64\"...")
            progress_bar.progress((i + 1) / len(choke_sizes))
            
            # 1. Find Operating Rate for this Choke
            # We assume Natural Flow intersection logic
            q_guess = 3000
            
            # Simple iteration to find rate where WHP_available == WHP_required
            # (Simplified solver for speed)
            for _ in range(5):
                bhp = f_ipr(q_guess)
                whp_avail = VLPModels.calc_whp_from_bhp(bhp, vlp['depth'], vlp['id'], q_guess, 
                                                        pvt['gor'], vlp['wc'], pvt['api'], pvt['yg'], pvt['t'])
                whp_req = ChokeModels.calc_pwh(q_guess, size, pvt['gor'], "Ros")
                
                diff = whp_avail - whp_req
                if abs(diff) < 10: break
                
                # Adjust guess
                if diff > 0: q_guess += 100 # Can flow more
                else: q_guess -= 100 # Choke restricted
            
            final_q = max(0, q_guess)
            
            # 2. Check Constraints at this Rate
            # Velocity Check
            _, _, _, v_act, v_eros_limit = VLPModels.gradient_traverse(
                 vlp['depth'], vlp['id'], whp_avail, final_q, 
                 pvt['gor'], pvt['api'], pvt['yg'], vlp['wc'], pvt['t']
            )
            max_v_in_tubing = max(v_act) if v_act else 0
            
            # Store Result
            status = "✅ Safe"
            if max_v_in_tubing > max_erosional_vel: status = "❌ Erosion Risk"
            if whp_avail < min_whp: status = "❌ Low Pressure"
            if final_q < 10: status = "❌ No Flow"
            
            results.append({
                "Choke (/64)": size,
                "Rate (STB/D)": final_q,
                "WHP (psi)": whp_avail,
                "Velocity (ft/s)": max_v_in_tubing,
                "Status": status
            })
            
        # 4. Display Results
        df_opt = pd.DataFrame(results)
        
        # Find Best Safe Option
        safe_opts = df_opt[df_opt['Status'] == "✅ Safe"]
        if not safe_opts.empty:
            best_opt = safe_opts.loc[safe_opts['Rate (STB/D)'].idxmax()]
            
            st.success(f"""
            ### 🏆 Optimal Configuration Found:
            * **Choke Size:** {best_opt['Choke (/64)']}/64"
            * **Max Safe Production:** {best_opt['Rate (STB/D)']:.0f} STB/D
            * **WHP:** {best_opt['WHP (psi)']:.0f} psi
            """)
        else:
            st.error("No configuration met all safety constraints. Consider increasing tubing size.")
            
        # Table with Color Formatting
        st.write("### Simulation Results")
        
        def color_status(val):
            color = 'green' if val == "✅ Safe" else 'red'
            return f'color: {color}'
            
        st.dataframe(df_opt.style.map(color_status, subset=['Status']), use_container_width=True)
        
        # Plot
        fig, ax = plt.subplots()
        ax.plot(df_opt['Choke (/64)'], df_opt['Rate (STB/D)'], 'o-')
        ax.axhline(y=safe_opts['Rate (STB/D)'].max() if not safe_opts.empty else 0, color='g', linestyle='--', label="Max Safe Rate")
        ax.set_xlabel("Choke Size (/64\")")
        ax.set_ylabel("Production Rate (STB/D)")
        ax.set_title("Choke Performance Optimization")
        ax.legend()
        ax.grid(True)
        st.pyplot(fig)
# --- MODULE 8: ARTIFICIAL LIFT & DIAGNOSTICS ---
elif module == "8. Artificial Lift & Diagnostics":
    st.title("Step 8: Artificial Lift & Diagnostics Studio")
    
    if not st.session_state.ipr or not st.session_state.vlp:
        st.error("⚠️ Complete Steps 2 & 3 first.")
        st.stop()
        
    # --- TABS ---
    tab_lift, tab_doc, tab_eco = st.tabs(["🏗️ Artificial Lift Simulator", "🩺 Well Doctor (Diagnostics)", "💰 Economic Optimizer"])
    
    # === TAB 1: ARTIFICIAL LIFT SIMULATOR ===
    with tab_lift:
        st.subheader("Simulate Lift Methods")
        st.write("Compare Natural Flow vs. Gas Lift vs. Pumps (ESP/SRP).")
        
        # 1. Select Lift Method
        lift_type = st.radio("Select Lift Method:", ["Natural Flow", "Gas Lift", "ESP / Downhole Pump"], horizontal=True)
        
        lift_val = 0
        if lift_type == "Gas Lift":
            lift_val = st.slider("Gas Injection Rate (Mscf/d)", 0, 5000, 1000)
            # Convert Mscf/d to approx Injection GLR for simple model
            # Assuming average rate of 2000 bbl/d for GLR estimation
            est_inj_glr = (lift_val * 1000) / 2000 
            st.caption(f"Estimated Injection GLR: +{est_inj_glr:.0f} scf/stb")
            param = est_inj_glr
        elif lift_type == "ESP / Downhole Pump":
            lift_val = st.slider("Pump Delta P Boost (psi)", 0, 3000, 500)
            st.caption("Pressure added by pump to overcome head.")
            param = lift_val
        else:
            param = 0

        # 2. Run Simulation
        q_scan = np.linspace(100, 6000, 60)
        
        # Base IPR
        ipr = st.session_state.ipr
        df_ipr = pd.DataFrame({'q': ipr['q'], 'p': ipr['p']}).sort_values('q')
        f_ipr = interp1d(df_ipr['q'], df_ipr['p'], fill_value="extrapolate")
        
        # Calculate VLP with Lift
        vlp = st.session_state.vlp
        pvt = st.session_state.pvt
        
        vlp_lift = []
        for q in q_scan:
            p_bh = VLPModels.calc_bh_pressure(
                q, vlp['thp'], vlp['depth'], vlp['id'], 
                pvt['gor'], vlp['wc'], pvt['api'], pvt['yg'], pvt['t'],
                lift_method=("Gas Lift" if lift_type=="Gas Lift" else ("ESP/Pump" if lift_type=="ESP / Downhole Pump" else "None")),
                lift_param=param
            )
            vlp_lift.append(p_bh)
            
        # 3. Find Intersection
        f_vlp_lift = interp1d(q_scan, vlp_lift, fill_value="extrapolate")
        diff = np.abs(f_ipr(q_scan) - vlp_lift)
        idx = np.argmin(diff)
        op_rate, op_press = q_scan[idx], f_ipr(q_scan)[idx]
        
        # 4. Plot
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(df_ipr['q'], df_ipr['p'], 'g-', linewidth=3, label="IPR")
        
        # Plot Natural VLP for reference
        if lift_type != "Natural Flow":
            vlp_nat = [VLPModels.calc_bh_pressure(q, vlp['thp'], vlp['depth'], vlp['id'], pvt['gor'], vlp['wc'], pvt['api'], pvt['yg'], pvt['t']) for q in q_scan]
            ax.plot(q_scan, vlp_nat, 'k--', alpha=0.5, label="Natural VLP")
            
        ax.plot(q_scan, vlp_lift, 'r-', linewidth=3, label=f"VLP ({lift_type})")
        
        if diff[idx] < 200:
            ax.plot(op_rate, op_press, 'bo', markersize=15)
            st.metric("Production Rate", f"{op_rate:.0f} STB/D", delta=None if lift_type=="Natural Flow" else f"Lift Effect: {lift_type}")
        else:
            st.error("No Intersection (Well Dead)")
            
        ax.set_xlabel("Rate (STB/D)"); ax.set_ylabel("Pressure (psi)")
        ax.legend(); ax.grid(True)
        st.pyplot(fig)

    # === TAB 2: WELL DOCTOR ===
    with tab_doc:
        st.subheader("🩺 Diagnostic Analysis")
        
        if diff[idx] < 200:
            col_d1, col_d2 = st.columns(2)
            
            # 1. Liquid Loading Check (Turner)
            # Turner Velocity (simplified): v_crit ~ 10-15 ft/s usually safe
            # We calculate actual velocity at wellhead
            _, _, _, v_act_profile, _ = VLPModels.gradient_traverse(
                vlp['depth'], vlp['id'], vlp['thp'], op_rate, 
                pvt['gor'], pvt['api'], pvt['yg'], vlp['wc'], pvt['t']
            )
            v_surface = v_act_profile[0] # Velocity at surface
            
            with col_d1:
                st.write("**1. Liquid Loading Check**")
                st.metric("Gas Velocity at Surface", f"{v_surface:.1f} ft/s")
                if v_surface < 10:
                    st.error("⚠️ CRITICAL: Velocity < 10 ft/s. High risk of liquid loading/well dying.")
                elif v_surface < 15:
                    st.warning("⚠️ WARNING: Velocity low. Monitor closely.")
                else:
                    st.success("✅ STABLE: Velocity sufficient to lift liquids.")
            
            # 2. Stability Check (Slope Analysis)
            # If VLP slope is negative at intersection, flow is unstable
            slope_vlp = (vlp_lift[idx+1] - vlp_lift[idx-1]) / (q_scan[idx+1] - q_scan[idx-1])
            
            with col_d2:
                st.write("**2. Flow Stability**")
                if slope_vlp < 0:
                    st.error(f"⚠️ UNSTABLE: VLP Slope is Negative ({slope_vlp:.2f}). Heading/Slugging likely.")
                else:
                    st.success(f"✅ STABLE: VLP Slope is Positive ({slope_vlp:.2f}).")
            
            # 3. Reservoir Health
            st.write("**3. Reservoir Pressure Drawdown**")
            dd_percent = (ipr['pr'] - op_press) / ipr['pr'] * 100
            st.progress(min(1.0, dd_percent/100))
            st.caption(f"Drawdown: {dd_percent:.1f}% of Reservoir Pressure")
            
        else:
            st.warning("Establish an operating point in the Simulator tab first.")

    # === TAB 3: ECONOMIC OPTIMIZER ===
    with tab_eco:
        st.subheader("💰 Economic Calculator")
        
        ce1, ce2, ce3 = st.columns(3)
        oil_price = ce1.number_input("Oil Price ($/bbl)", value=75.0)
        water_cost = ce2.number_input("Water Disposal ($/bbl)", value=2.0)
        lift_cost = ce3.number_input("Lift Cost ($/day)", value=100.0) # Electricity/Gas
        
        if diff[idx] < 200:
            # Calculate Revenue
            q_oil = op_rate * (1 - vlp['wc'])
            q_water = op_rate * vlp['wc']
            
            rev_oil = q_oil * oil_price
            cost_water = q_water * water_cost
            
            daily_profit = rev_oil - cost_water - lift_cost
            
            st.divider()
            c_res1, c_res2, c_res3 = st.columns(3)
            c_res1.metric("Daily Revenue", f"${rev_oil:,.0f}")
            c_res2.metric("Daily Costs", f"${(cost_water + lift_cost):,.0f}")
            c_res3.metric("Net Profit", f"${daily_profit:,.0f}", delta_color="normal")
            
            # Profit Breakdown Chart
            fig_eco, ax_eco = plt.subplots(figsize=(6, 2))
            ax_eco.barh(['Profit', 'Costs'], [daily_profit, cost_water+lift_cost], color=['green', 'red'])
            st.pyplot(fig_eco)
        else:
            st.info("Run simulation to see economics.")