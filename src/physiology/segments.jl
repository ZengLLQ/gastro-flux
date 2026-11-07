"""
    Segments: Human Anatomical Gastrointestinal Compartments (ACAT Framework)
"""

"""
    create_human_physiology(; state::Symbol = :Fasted)

Constructs human anatomical and physiological state with 9 compartmental segments:
1. Stomach
2. Duodenum
3. Jejunum 1
4. Jejunum 2
5. Ileum 1
6. Ileum 2
7. Ileum 3
8. Caecum
9. Colon

Supports `:Fasted` and `:Fed` prandial states with corresponding pH, fluid volume,
transit time, and bile salt micellar concentrations.
"""
function create_human_physiology(; state::Symbol = :Fasted)
    is_fasted = (state == :Fasted)

    if is_fasted
        # Fasted state parameters
        stomach_vol = 50.0
        stomach_transit = 0.25 # h (~15 min gastric half-emptying)
        stomach_pH = 1.7
        stomach_bile = 0.05    # mM

        duo_vol = 45.0
        duo_transit = 0.26
        duo_pH = 6.0
        duo_bile = 3.5

        jej1_vol = 105.0
        jej1_transit = 0.55
        jej1_pH = 6.5
        jej1_bile = 2.8

        jej2_vol = 105.0
        jej2_transit = 0.55
        jej2_pH = 6.8
        jej2_bile = 2.0

        ile1_vol = 75.0
        ile1_transit = 0.50
        ile1_pH = 7.2
        ile1_bile = 1.5

        ile2_vol = 75.0
        ile2_transit = 0.50
        ile2_pH = 7.5
        ile2_bile = 1.2

        ile3_vol = 75.0
        ile3_transit = 0.50
        ile3_pH = 7.5
        ile3_bile = 1.0

        caecum_vol = 45.0
        caecum_transit = 4.5
        caecum_pH = 6.4
        caecum_bile = 0.5

        colon_vol = 150.0
        colon_transit = 18.0
        colon_pH = 6.8
        colon_bile = 0.1

        t_half_stomach = 0.25
    else
        # Fed state parameters
        stomach_vol = 250.0
        stomach_transit = 1.5 # Delayed meal gastric emptying
        stomach_pH = 4.5      # Food buffer
        stomach_bile = 1.5

        duo_vol = 65.0
        duo_transit = 0.35
        duo_pH = 5.4
        duo_bile = 12.0       # Postprandial gallbladder bile release

        jej1_vol = 140.0
        jej1_transit = 0.65
        jej1_pH = 6.2
        jej1_bile = 9.5

        jej2_vol = 140.0
        jej2_transit = 0.65
        jej2_pH = 6.6
        jej2_bile = 6.0

        ile1_vol = 90.0
        ile1_transit = 0.55
        ile1_pH = 7.0
        ile1_bile = 4.0

        ile2_vol = 90.0
        ile2_transit = 0.55
        ile2_pH = 7.4
        ile2_bile = 2.5

        ile3_vol = 90.0
        ile3_transit = 0.55
        ile3_pH = 7.5
        ile3_bile = 1.5

        caecum_vol = 60.0
        caecum_transit = 5.0
        caecum_pH = 6.5
        caecum_bile = 0.8

        colon_vol = 200.0
        colon_transit = 20.0
        colon_pH = 6.8
        colon_bile = 0.2

        t_half_stomach = 1.5
    end

    # Construct 9 anatomical segments
    segments = [
        GISegment("Stomach",   stomach_vol, stomach_transit, stomach_pH, stomach_bile, 20.0,  3.5,  1.0, 0.00, 0.00),
        GISegment("Duodenum",  duo_vol,     duo_transit,     duo_pH,     duo_bile,     25.0,  1.5,  4.0, 0.35, 0.20),
        GISegment("Jejunum 1", jej1_vol,    jej1_transit,    jej1_pH,    jej1_bile,    100.0, 1.25, 4.0, 1.00, 0.50),
        GISegment("Jejunum 2", jej2_vol,    jej2_transit,    jej2_pH,    jej2_bile,    100.0, 1.25, 3.5, 0.80, 0.70),
        GISegment("Ileum 1",   ile1_vol,    ile1_transit,    ile1_pH,    ile1_bile,    100.0, 1.0,  2.5, 0.50, 0.90),
        GISegment("Ileum 2",   ile2_vol,    ile2_transit,    ile2_pH,    ile2_bile,    100.0, 1.0,  2.0, 0.30, 1.00),
        GISegment("Ileum 3",   ile3_vol,    ile3_transit,    ile3_pH,    ile3_bile,    100.0, 1.0,  1.5, 0.10, 1.00),
        GISegment("Caecum",    caecum_vol,  caecum_transit,  caecum_pH,  caecum_bile,  15.0,  3.5,  1.0, 0.00, 0.80),
        GISegment("Colon",     colon_vol,   colon_transit,   colon_pH,   colon_bile,   120.0, 2.5,  1.0, 0.00, 0.60)
    ]

    # Systemic disposition parameters: Q_H = 90 L/h, Q_portal = 72 L/h
    # 3-compartment body: V_c = 15 L, V_p1 = 25 L, V_p2 = 40 L
    return PhysiologyState(
        is_fasted,
        t_half_stomach,
        segments,
        90.0,   # liver_blood_flow_L_h
        72.0,   # portal_blood_flow_L_h
        15.0,   # systemic_v_central_L
        25.0,   # systemic_v_periph1_L
        40.0,   # systemic_v_periph2_L
        0.35,   # systemic_k12_h
        0.22,   # systemic_k21_h
        0.06,   # systemic_k13_h
        0.02    # systemic_k31_h
    )
end
