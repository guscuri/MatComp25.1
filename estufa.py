import numpy as np
import requests
import matplotlib.pyplot as plt
#esturuta com 40x40x80cm

# --- Parâmetros do sistema ---
v_est = 0.128 # volume da estufa (m³)
k_est = 0.036 # coeficiente de troca térmica (W/m²·K)
p_ar = 1.225 # densidade do ar (kg/m³)
A_est = 1.44 # área da estufa (m²)
d = 0.0001
v_dot = 0.0001 # vazão de ar do ventilador (m³/s)

# Inércia térmica do ar
m_ar = p_ar * v_est    # massa de ar (kg)
c_p = 1005          # calor específico do ar (J/kg·K)
C_T = m_ar * c_p    # capacidade térmica (J/K)

# Coeficiente global de troca térmica (WA)
U = k_est/d           # W/K*m²
UA = U * A_est       # W/K

# Potência do LED
Q_led = 80    # W

# Renovação de ar
m_dot = p_ar * v_dot
m_dot_cp = m_dot * c_p    # vazão mássica vezes c_p (W/K)
# Umidade relativa (Fazer teste na estufa para descorbrir valores reais)
k_evap = 0.0005     # coef. evap. natural (kg/s·m³) 
k_hum = 0.0002      # ganho umidificador (kg/s/unit)
k_desum = 0.0002    # ganho desumidificador (kg/s/unit)


# Condições externas e de entrada (podem ser funções de t)
# Get real-time temperature for Copacabana, Rio de Janeiro

API_KEY = "512bd22151f6d8c4f75fe2ae248168c1"  # Replace with your actual API key
LOCATION = "Copacabana,Rio de Janeiro,BR"

def get_copacabana(x):
    try:
        url = f"http://api.openweathermap.org/data/2.5/weather?q={LOCATION}&appid={API_KEY}&units=metric"
        response = requests.get(url)
        data = response.json()
        if x == "temp":
            return data['main']['temp']
        elif x == "humi":
            return data['main']['humidity']
        else:
            raise ValueError("Invalid parameter. Use 'temp' or 'humi'.")
    except Exception as e:
        temp = data['main']['temp']
        humidity = data['main']['humidity']  # Get humidity percentage
        #print(f"Current humidity in Copacabana: {humidity}%")
        return temp, humidity
    except Exception as e:
        print(f"Error fetching temperature: {e}")
        return 25.0  # Default fallback temperature

# Get current temperature from API
T_amb = get_copacabana("temp")  # temperatura ambiente externa (°C)
T_ent = T_amb           # temperatura do ar de entrada (°C)
T_int = 23.0         # temperatura interna inicial (°C)
T_obj_min, T_obj_max = 18.0, 26.0      # intervalo de temperatura interna desejada (°C)
U_amb = get_copacabana("humi")  # umidade relativa ambiente externa (%)
#print(f"Current temperature in Copacabana: {T_amb}°C")

# --- Definição da EDO ---
# --- Temperatura da estufa ---
def dT_dt(T, t, T_amb, T_ent):
    """
    Retorna dT/dt para o modelo de temperatura da estufa.
    """
    term_conducao = UA * (T_amb - T)
    term_led     = Q_led
    term_vent    = m_dot_cp * (T_ent - 2*T + T_amb)
    return (term_conducao + term_led + term_vent) / C_T
# --- Umidade Relativa da estufa ---

# Saturação (Clausius–Clapeyron)
def H_sat(T):
    # T em °C, retornar H_sat em kg/m³
    e_sat = 0.611 * np.exp(17.27 * T / (T + 237.3))  # kPa
    return 216.7 * e_sat / (T + 273.15)

def dH_dt(H, T, u_hum, u_desum):
    term_evap  = k_evap * (H_sat(T) - H)
    term_hum   = k_hum * u_hum
    term_desum = k_desum * u_desum
    return (term_evap + term_hum - term_desum) / v_est

# --- Métodos RK4 ---
# ---Temperatura da estufa ---
def rk4_step_t(fun, T, t, dt, *args):
    k1 = fun(T, t, *args)
    k2 = fun(T + 0.5*dt*k1, t + 0.5*dt, *args)
    k3 = fun(T + 0.5*dt*k2, t + 0.5*dt, *args)
    k4 = fun(T + dt*k3,   t + dt,   *args)
    return T + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

# --- Umidade Relativa da estufa ---
# RK4 genérico para umidade
def rk4_step_h(fun, H, T, dt, u_hum, u_desum):
    k1 = fun(H, T, u_hum, u_desum)
    k2 = fun(H + 0.5*dt*k1, T, u_hum, u_desum)
    k3 = fun(H + 0.5*dt*k2, T, u_hum, u_desum)
    k4 = fun(H +     dt*k3, T, u_hum, u_desum)
    return H + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

"""# --- Simulação aberta (sem controle) ---
t_end = 3600       # 1 hora em segundos
dt = 1.0           # passo de integração (s)
nt = int(t_end/dt)

# Vetores de tempo e temperatura
ts = np.linspace(0, t_end, nt+1)
Ts = np.zeros(nt+1)
Ts[0] = T_amb      # condição inicial: igual ao ambiente

# Loop de simulação
t = 0.0
for i in range(nt):
    Ts[i+1] = rk4_step_t(dT_dt, Ts[i], t, dt, T_amb, T_ent)
    t += dt

# --- Plot dos resultados ---
plt.figure()
plt.plot(ts/60, Ts, label='Temperatura interna (°C)')
plt.axhline(T_amb, color='gray', linestyle='--', label='T ambiente')
plt.xlabel('Tempo (min)')
plt.ylabel('T (°C)')
plt.title('Simulação de Temperatura da Estufa (aberta)')
plt.legend()
plt.grid(True)
plt.show()
"""
