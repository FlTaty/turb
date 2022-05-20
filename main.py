from iapws import IAPWS97 as WSP
from iapws import IAPWS97 as IAPWS97 
import math as M
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
import pandas as pd
from sympy import *


st.write('Луговых Т. ФПэ-01-19')
st.write("Github: " + "https://github.com/tederix/PGT")

page = st.sidebar.selectbox(
    "Выберите задание",
    ("Задание 1", "Задание 2", "Задание 3"))

if page == "Задание 1":
    st.write("""# """)
    st.write("""Построить процесс расширения пара в турбине. Определение расходов пара на входе в турбину (G0) и в конденсатор (Gк). Получить зависимость КПД ПТУ от параметра заданного в таблице.""")
    st.write("""# """)
    st.write(" *Исходные данные:* ")
    Ne = st.number_input('Введите мощность Nэ, МВт', value = 830)*10**6
    t0 = st.number_input('Введите температуру T0, °C', value = 560)
    T0 = t0+273.15
    ppp = st.number_input('Введите давление P0, МПа', value = 3.62)*10**6
    tpp = st.number_input('Введите температуру Tпп, °C', value = 565)
    Tpp = tpp+273.15
    pk = st.number_input('Введите давление Pk, кПа', value = 3.5)*10**3
    tpv = st.number_input('Введите температуру Tпв, °C', value = 269)
    Tpv = tpv+273.15

    age = st.slider('Укажите максимальную границу P0', min_value = 20.0, max_value = 26.0, step = 1.0)
    age = age + 0.01
    P_0 = list(np.arange(20, age, 1))
    p0 = [p*1e6 for p in P_0]
    p_0_min = float(p0[0])
    p_0_max = float(p0[-1])

    z = 9

    st.write("""# """)
    st.write(" *Дано:* ")
    st.write(""" t0 = """ + str(t0) + """ C""")
    st.write(""" P0 = """ + str(p_0_min*10**(-6)) + " - " + str('{:.2}'.format(p_0_max*10**(-6))) + """ МПа""")
    st.write(""" tпп = """ + str(tpp) + """ C """)
    st.write(""" Pпп = """ + str(ppp*10**(-6)) + """ МПа""")
    st.write(""" Pк = """ + str(pk*10**(-3)) + """ кПа """)
    st.write(""" tпв = """ + str(tpv) + """ C """)
    st.write(""" Nэ = """ + str(Ne*10**(-6)) + """ МВт """)
    st.write(""" Z = """ + str(z) + """ шт """)

    st.write("""# """)
    st.write(" *Решение:* ")
    
    def Calculate_eta_G0_Gk(N_e, p_0, T_0, p_pp, T_pp, p_k, T_pv):

        point_0 = IAPWS97(P=p_0*10**(-6), T=T_0)
        s_0 = point_0.s
        h_0 = point_0.h
        v_0 = point_0.v
        p_0_ = p_0-0.05*p_0
        point_p_0_ = IAPWS97(P=p_0_*10**(-6), h=h_0)
        t_0_ = point_p_0_.T-273.15
        s_0_ = point_p_0_.s
        v_0_ = point_p_0_.v

        p_1t = p_pp+0.1*p_pp
        point_1t = IAPWS97(P=p_1t*10**(-6), s=s_0)
        t_1t = point_1t.T-273.15
        h_1t = point_1t.h
        v_1t = point_1t.v

        point_pp = IAPWS97 (P=p_pp*10**(-6), T=T_pp)
        h_pp = point_pp.h
        s_pp = point_pp.s
        v_pp = point_pp.v

        H_0 = h_0-h_1t
        eta_oi = 0.85
        H_i_cvd = H_0*eta_oi

        h_1 = h_0 - H_i_cvd
        point_1 = IAPWS97(P = p_1t*10**(-6),h = h_1)
        s_1 = point_1.s
        T_1 = point_1.T
        v_1 = point_1.v
        p_pp_ = p_pp - 0.03*p_pp
        point_pp_ = IAPWS97(P=p_pp_*10**(-6),h = h_pp)
        s_pp_ = point_pp_.s
        v_pp_ = point_pp_.v
        point_kt = IAPWS97(P = p_k*10**(-6),s = s_pp)
        T_kt = point_kt.T
        h_kt = point_kt.h
        v_kt = point_kt.v
        s_kt = s_pp
        H_0_csdcnd = h_pp-h_kt
        eta_oi = 0.85
        H_i_csdcnd = H_0_csdcnd*eta_oi
        h_k = h_pp - H_i_csdcnd
        point_k = IAPWS97(P = p_k*10**(-6), h = h_k)
        T_k = point_k.T
        s_k = point_k.s
        v_k = point_k.v
        point_k_v = IAPWS97(P = p_k*10**(-6),x=0)
        h_k_v = point_k_v.h
        s_k_v = point_k_v.s
        eta_oiI = (h_1-h_0)/(h_1t-h_0)
        p_pv = 1.4*p_0
        point_pv = IAPWS97(P = p_pv*10**(-6),T=T_pv)
        h_pv = point_pv.h
        s_pv = point_pv.s
        ksi_pp_oo = 1-(1-(T_k*(s_pp-s_k_v))/((h_0-h_1t)+(h_pp-h_k_v)))/(1-(T_k*(s_pp-s_pv))/((h_0-h_1t)+(h_pp-h_pv)))
        #T_0_= IAPWS97(P = p_pv*10**(-6),x = 0).T
        T_0_ = 374.2+273.15
        T_ = (point_pv.T - point_k.T) / (T_0_ - point_k.T)
        if T_ <= 0.636364:
            ksi1 = -1.53*T_**2+2.1894*T_+0.0048
        elif 0.636364<T_<=0.736364:
            ksi1 = -1.3855*T_**2+2.0774*T_+0.0321
        elif 0.736364<T_<=0.863636:
            ksi1 = -2.6535*T_**2+4.2556*T_-0.8569

        if T_ <= 0.636364:
            ksi2 = -1.53*T_**2+2.1894*T_+0.0048
        elif 0.636364<T_<=0.736364:
            ksi2 = -2.5821*T_**2+3.689*T_-0.4825
        elif 0.718182<T_<=0.827273:
            ksi2 = -1.3855*T_**2+2.0774*T_+0.0321
        elif 0.736364<T_<=0.863636:
            ksi2 = -2.6535*T_**2+4.2556*T_-0.8569

        ksi = (ksi1+ksi2)/2
        ksi_r_pp = ksi*ksi_pp_oo
        eta_ir = (H_i_cvd+H_i_csdcnd)/(H_i_cvd+(h_pp-h_k_v))*1/(1-ksi_r_pp)
        H_i = eta_ir*((h_0-h_pv)+(h_pp-h_1))
        eta_m = 0.994
        eta_eg = 0.99
        G_0 = N_e/(H_i*eta_m*eta_eg*(10**3))
        G_k = N_e/((h_k-h_k_v)*eta_m*eta_eg*(10**3))*(1/eta_ir-1)

        return eta_ir, G_0, G_k

    eta, G0, Gk =[], [], []
    for p in p0:
        eta_ = Calculate_eta_G0_Gk(N_e = Ne, p_0 = p, T_0 = T0, p_pp = ppp, T_pp = Tpp, p_k = pk, T_pv = Tpv)
        eta.append(eta_[0])
        G0.append(eta_[1])
        Gk.append(eta_[2])

    max: float = eta[0]
    pos = 0
    for i in range(len(eta)):
        if eta[i] > max: max = eta[i]; pos = i

    delta_p_0 = 0.05*p0[pos]
    delta_p_pp = 0.08*ppp
    delta_p = 0.03*ppp

    p0_f = [float(x) * 10**(-6) for x in p0]
    eta_f = [float(x) * 100 for x in eta]

    st.write(""" Максимальное КПД = """ + str('{:.4}'.format(float(eta_f[pos]))) + """ %""")
    st.write(""" Расход пара на входе в турбину (G0) при макс. КПД = """ + str('{:.5}'.format(float(G0[pos]))) + """ кг/с""")
    st.write(""" Расход пара на входе в конденсатор (Gк) при макс. КПД = """ + str('{:.5}'.format(float(Gk[pos]))) + """ кг/с""")
    st.write("""# """)
    st.write(" Табл. Зависимость КПД от Pпп  ")

    p0_eta=pd.DataFrame({"p0, МПа": (p0_f),
                       "eta, %": (eta_f),
                       "G_0, кг/с": (G0),
                       "G_k, кг/с": (Gk)
                       })
    st.dataframe(p0_eta)

    st.write("""# """)

    p0__eta = plt.figure()

    plt.plot(p0_f, eta_f)
    plt.plot(p0_f, eta_f, 'ro')
    plt.title("Зависимость КПД от начального давления пара")
    plt.xlabel("P0, MПа")
    plt.ylabel("КПД, %")
    plt.grid()

    st.pyplot(p0__eta)

    st.title(""" """)

    p_0_max_max = p0[pos]

    fighs = plt.figure()

    point_0 = IAPWS97(P=p_0_max*1e-6, T=T0)
    p_0_d = p_0_max - delta_p_0
    point_0_d = IAPWS97(P=p_0_d*1e-6, h=point_0.h)
    p_1t = ppp + delta_p_pp
    point_1t = IAPWS97(P=p_1t*10**(-6), s=point_0.s)
    H_01 = point_0.h - point_1t.h
    kpd_oi = 0.85
    H_i_cvd = H_01 * kpd_oi
    h_1 = point_0.h - H_i_cvd
    point_1 = IAPWS97(P=p_1t*1e-6, h=h_1)
    point_pp = IAPWS97(P=ppp*1e-6, T=Tpp)
    p_pp_d = ppp - delta_p_pp
    point_pp_d = IAPWS97(P=p_pp_d*1e-6, h=point_pp.h)
    point_kt = IAPWS97(P=pk*1e-6, s=point_pp.s)
    H_02 = point_pp.h - point_kt.h
    kpd_oi = 0.85
    H_i_csd_cnd = H_02 * kpd_oi
    h_k = point_pp.h - H_i_csd_cnd
    point_k = IAPWS97(P=pk*1e-6, h=h_k)

    s_0 = [point_0.s-0.05,point_0.s,point_0.s+0.05]
    h_0 = [IAPWS97(P = p_0_max*1e-6,s = s_).h for s_ in s_0]
    s_1 = [point_0.s-0.05,point_0.s,point_0.s+0.18]
    h_1 = [IAPWS97(P=p_1t*1e-6, s = s_).h for s_ in s_1]
    s_0_d = [point_0_d.s-0.05, point_0_d.s, point_0_d.s+0.05]
    h_0_d = h_0
    s_pp = [point_pp.s-0.05,point_pp.s,point_pp.s+0.05]
    h_pp = [IAPWS97(P=ppp*1e-6, s=s_).h for s_ in s_pp]
    s_k = [point_pp.s-0.05,point_pp.s,point_pp.s+0.8]
    h_k = [IAPWS97(P=pk*1e-6, s=s_).h for s_ in s_k]
    s_pp_d = [point_pp_d.s-0.05,point_pp_d.s,point_pp_d.s+0.05]
    h_pp_d = h_pp

    plt.plot([point_0.s,point_0.s,point_0_d.s,point_1.s],[point_1t.h,point_0.h,point_0.h,point_1.h],'-or')
    plt.plot([point_pp.s,point_pp.s,point_pp_d.s,point_k.s],[point_kt.h,point_pp.h,point_pp.h,point_k.h],'-or')
    plt.plot(s_0,h_0)
    plt.plot(s_1,h_1)
    plt.plot(s_0_d,h_0_d)
    plt.plot(s_pp,h_pp)
    plt.plot(s_k,h_k)
    plt.plot(s_pp_d,h_pp_d)

    for x, y, ind in zip([point_pp.s, point_k.s], [point_pp.h, point_k.h], ['{пп}', '{к}']):
        plt.text(x-0.45, y+40, '$h_' + ind + ' = %.2f $'%y)
    for x, y, ind in zip([point_kt.s, point_pp_d.s], [point_kt.h, point_pp_d.h], ['{кт}', '{ппд}']):
        plt.text(x+0.03, y+40, '$h_' + ind + ' = %.2f $'%y)
    for x, y, ind in zip ([point_0.s, point_1.s], [point_0.h, point_1.h], ['{0}', '{1}']):
        plt.text(x-0.01, y+120, '$h_' + ind + ' = %.2f $'%y)
    for x, y, ind in zip([point_1t.s, point_0_d.s], [point_1t.h, point_0_d.h], ['{1т}', '{0д}']):
        plt.text(x+0.03, y-60, '$h_' + ind + ' = %.2f $'%y)

        plt.title("h - s диаграмма")
        plt.xlabel("s, кДж/(кг*С)")
        plt.ylabel("h, кДж/кг")
        plt.grid(True)


    st.pyplot(fighs)
    
