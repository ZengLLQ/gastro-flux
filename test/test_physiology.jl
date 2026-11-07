using Test
using GastroFlux

@testset "Gastrointestinal Physiology & Transporters" begin
    phys_fasted = create_human_physiology(state=:Fasted)
    phys_fed    = create_human_physiology(state=:Fed)

    # 1. Segment anatomy count and ordering
    @test length(phys_fasted.segments) == 9
    @test length(phys_fed.segments) == 9
    @test phys_fasted.segments[1].name == "Stomach"
    @test phys_fasted.segments[2].name == "Duodenum"
    @test phys_fasted.segments[9].name == "Colon"

    # 2. Fasted vs Fed state physiological differences
    # Stomach volume and pH
    @test phys_fasted.segments[1].volume_mL < phys_fed.segments[1].volume_mL
    @test phys_fasted.segments[1].pH < phys_fed.segments[1].pH
    @test phys_fasted.gastric_half_emptying_h < phys_fed.gastric_half_emptying_h

    # Gallbladder bile release: duodenum and jejunum bile salts are much higher in fed state
    @test phys_fed.segments[2].bile_salt_mM > phys_fasted.segments[2].bile_salt_mM
    @test phys_fed.segments[3].bile_salt_mM > phys_fasted.segments[3].bile_salt_mM

    # 3. Active P-gp efflux transporter saturation
    drug_pgp = DrugProperties(
        name="PgpSubstrate",
        vmax_pgp=5.0, # mg/h
        km_pgp=0.01   # mg/mL
    )
    seg_jej = phys_fasted.segments[3]
    vol_ent = 5.0 # mL

    # Low concentration: rate proportional to concentration
    rate_low = efflux_rate(drug_pgp, seg_jej, 0.001, vol_ent)
    # High concentration: approaches Vmax * pgp_rel
    rate_high = efflux_rate(drug_pgp, seg_jej, 10.0, vol_ent)

    @test rate_low > 0.0
    @test rate_high > rate_low
    max_theoretical = drug_pgp.vmax_pgp * seg_jej.pgp_relative
    @test rate_high < max_theoretical
    @test isapprox(rate_high, max_theoretical, rtol=0.05)

    # 4. Enterocyte CYP3A4 metabolism saturation
    drug_cyp = DrugProperties(
        name="CypSubstrate",
        vmax_cyp3a4=12.0, # mg/h
        km_cyp3a4=0.02    # mg/mL
    )
    cyp_rate_low = gut_metabolism_rate(drug_cyp, seg_jej, 0.002, vol_ent)
    cyp_rate_high = gut_metabolism_rate(drug_cyp, seg_jej, 20.0, vol_ent)

    @test cyp_rate_low > 0.0
    @test cyp_rate_high > cyp_rate_low
    max_cyp_theoretical = drug_cyp.vmax_cyp3a4 * seg_jej.cyp3a4_relative
    @test isapprox(cyp_rate_high, max_cyp_theoretical, rtol=0.05)

    # 5. Hepatic disposition well-stirred model
    drug_high_cl = DrugProperties(hepatic_clearance_int=500.0, fraction_unbound_plasma=0.5)
    cl_h, e_h, f_h = hepatic_disposition(drug_high_cl, phys_fasted.liver_blood_flow_L_h)

    @test 0.0 <= e_h <= 1.0
    @test 0.0 <= f_h <= 1.0
    @test isapprox(e_h + f_h, 1.0, atol=1e-12)
    @test cl_h < phys_fasted.liver_blood_flow_L_h # Blood flow limited
end
