function resolution_modele(instance)
    # Pré-traitement
    choix_pretraitement = choix_binaire("\n --> Souhaitez-vous désactiver les positions non-couvrantes (o/n) ? ")
    if choix_pretraitement == "o"
        temps_pretraitement = @elapsed positions_non_couvrantes(instance)
    end

    # Pré-calculs
    temps_precalcul_unites = @elapsed E1 = calcul_positions_unites(instance)
    temps_precalcul_bases = @elapsed E2 = calcul_positions_bases(instance)
    temps_precalcul_bases = @elapsed E4 = calcul_positions_porte_bases(instance)
    # Pré-calculs essentiels (niveau 2 : messages)
    Ebis = Int[] # On cherche les positions actives (sous-ensemble de positions)
    for position in 1:instance.contexte.nE
        if instance.contexte.Ebool[position] == true # position active
            push!(Ebis, position)

        end
    end
    temps_precalcul_positions_relais = @elapsed E3 = calcul_positions_relais(instance, Ebis) # positions à portée d'une autre position pour un relais donné

    # Création du modèle
    println("\n========== Création du modèle ==========")
    choix = choix_binaire("\n --> Souhaitez-vous ajouter la contrainte additionnelle R1 (o/n) ? ")
    ajout_R1 = choix == "o" ? true : false
    choix = choix_binaire("\n --> Souhaitez-vous ajouter la contrainte additionnelle R2 (o/n) ? ")
    ajout_R2 = choix == "o" ? true : false
    relaxed=true
    temp_creation_modele = @elapsed relaxed_modele = creation_modele(instance, E1, E2, E4, ajout_R1, ajout_R2, relaxed)
    G,DG,Vertices, Edges, Arcs=graph_conctruction(instance, Ebis)
    relaxed=false
    temp_creation_modele = @elapsed modele = creation_modele(instance, E1, E2, E4, ajout_R1, ajout_R2, relaxed)


    println("\n - Nombre de variables : ", num_variables(modele))
    nb_contraintes = 0
    for (F, S) in list_of_constraint_types(modele)
        nb_contraintes += num_constraints(modele, F, S)
    end
    println("\n - Nombre de contraintes : ", nb_contraintes)
    choix = choix_binaire("\n --> Souhaitez-vous enregistrer le modèle dans un fichier (o/n) ? ")
    if choix == "o"
        nom_dossier = string(today())*"_"*Dates.format(now(), "HH:MM:SS")
        run(`mkdir modeles/$(nom_dossier)`)
        write_to_file(modele, "modeles/$(nom_dossier)/modele.dat", format=MOI.FileFormats.FORMAT_LP)
        println("\n /!\\ Le modèle est disponible dans le dossier \"modeles/$(nom_dossier)/\" /!\\ ")
    end
    println("\n========== Création du modèle ==========")

    # Choix du solveur
    println("\n Choix du solveur : ")
    println(" ------------------ ")
    println("    Non-commerciaux ")
    println("    1) GLPK (GNU Linear Programming Kit)")
    println("    2) CBC (COIN-OR Branch and Cut)")
    println("    3) SCIP (Solving Constraint Integer Programs)")
    println("    Commerciaux ")
    println("    4) Gurobi")
    println("    5) CPLEX")
    println("    6) Mosek")
    choix = choix_multiple("\n --> Votre choix : ", 6)
    if choix == 1
        set_optimizer(modele, GLPK.Optimizer)
        choix_verbosite = gestion_parametres_GLPK(modele)
    elseif choix == 2
        set_optimizer(modele, Cbc.Optimizer)
        choix_verbosite = gestion_parametres_CBC(modele)
    elseif choix == 3
         set_optimizer(modele, SCIP.Optimizer)
         choix_verbosite = gestion_parametres_SCIP(modele)
    elseif choix == 4
        print("\n - ")
        set_optimizer(modele, Gurobi.Optimizer)
        choix_verbosite = gestion_parametres_Gurobi(modele)
    elseif choix == 5
        set_optimizer(modele, CPLEX.Optimizer)
        choix_verbosite = gestion_parametres_CPLEX(modele)
    elseif choix == 6
        set_optimizer(modele, Mosek.Optimizer)
        choix_verbosite = gestion_parametres_Mosek(modele)
    end

    # Limite de temps
    choix_p1 = choix_binaire("\n --> Souhaitez-vous fixer une limite de temps (o/n) ? ")
    if choix_p1 == "o"
        print("\n --> Temps limite (en secondes) : ")
        #temps_limite = parse(Float64, readline()) * 10^3 # Secondes vers milli-secondes
        #set_optimizer_attribute(modele, "tm_lim", temps_limite)
        temps_limite = parse(Float64, readline())
        set_time_limit_sec(modele, temps_limite)
    end

    println("solving the relaxation")
    set_optimizer(relaxed_modele, CPLEX.Optimizer)
    # relax_integrality(relaxed_modele)
    optimize!(relaxed_modele)
    z2 = objective_value(relaxed_modele)
    lecture_resultats_modele(instance, E1, E2, relaxed_modele, Vertices, Arcs, Edges, Ebis, true)
    println("relxation solution: $z2")

    #
    # Résolution
    if choix_verbosite == "o"
        println("\n========== Détails de la résolution ==========\n")
        optimize!(modele)
        print("\n========== Détails de la résolution ==========")
    else
        print("\n - Résolution en cours...")
        optimize!(modele)
    end


    # Affichage des temps
    precision = 4 # Précision des résultats (nombre de décimales)
    temps_resolution = solve_time(modele)
    println("\n\n========== Affichage des temps ==========")
    if choix_pretraitement == "o"
        println("\n - Temps pré-traitement : $(trunc(temps_pretraitement, digits=precision)) seconde(s)")
    end
    println("\n - Temps pré-calcul n°1 (unités) : $(trunc(temps_precalcul_unites, digits=precision)) seconde(s)")
    println("\n - Temps pré-calcul n°2 (bases) : $(trunc(temps_precalcul_bases, digits=precision)) seconde(s)")
    println("\n - Temps création modèle : $(trunc(temp_creation_modele, digits=precision)) seconde(s)")
    println("\n - Temps résolution : $(trunc(temps_resolution, digits=precision)) seconde(s)")
    println("\n========== Affichage des temps ==========")

    if termination_status(modele) == MOI.INFEASIBLE
        println("\n - Échec : problème non réalisable")
        statut = 0
    elseif termination_status(modele) == MOI.OPTIMAL # Solution optimale
        println("\n - Solution optimale trouvée")
        statut = 1
    elseif termination_status(modele) == MOI.TIME_LIMIT && has_values(modele) # Temps limite atteint mais solution réalisable sous-optimale disponible
        println("\n - Solution réalisable (possiblement sous-optimale) trouvée dans le temps imparti")
        statut = 1
    elseif termination_status(modele) == MOI.TIME_LIMIT # Temps limite atteint et pas de solution
        println("\n - Échec : temps limite de résolution atteint (> $(temps_limite/(10^3)) seconde(s)) ")
        statut = 0
    end

    if statut == 1 # Solution disponible
        print("\n --> Appuyez sur la touche \"Entrée\" pour afficher les résultats")
        readline()
        Ebis = Int[] # On cherche les positions actives (sous-ensemble de positions)
        for position in 1:instance.contexte.nE
            if instance.contexte.Ebool[position] == true # position active
                push!(Ebis, position)

            end
        end
        G,DG,Vertices, Edges, Arcs=graph_conctruction(instance, Ebis)

        Ub=calcul_unites_couverte_bases(instance)
        U=[(u,t) for u in 1:instance.contexte.nU, t in 1:instance.contexte.nT if (u,t) ∉  Ub]

        solution = lecture_resultats_modele(instance, E1, E2, modele, Vertices, Arcs, Edges, Ebis, relaxed)
    elseif statut == 0 # Pas de solution
        solution = Solution(0, 0, Int[], Int[], Int[], Int[], Int[])
    end

    z1=solution.z
    if z1==0
        Gap=0
    else
        println("solopt: $z1")
        println("solFractionnaire: $z2")
        Gap= ((z2-z1)/z1)*100
        # Gap= ((z1-z2)/z1)*100
    end

     println("Gap: $Gap")

    #solve the relaxation
    return solution
end



function creation_modele(instance, E1, E2, E4, ajout_R1, ajout_R2, relaxed)
    # Initialisation du modèle
    modele = Model()
    Ebis = Int[] # On cherche les positions actives (sous-ensemble de positions)
    Ebisbar = Int[] # On cherche les positions actives (sous-ensemble de positions)

    for position in 1:instance.contexte.nE
        if instance.contexte.Ebool[position] == true # position active
            push!(Ebis, position)
        else
            push!(Ebisbar, position)

        end
    end


    E3 = calcul_positions_relais(instance, Ebis) # positions à portée d'une autre position pour un relais donné

    # Données
    nT = instance.contexte.nT
    nH = instance.contexte.nH
    nE = instance.contexte.nE
    nC = instance.contexte.nC
    nU = instance.contexte.nU
    nB = instance.contexte.nB
    nR = instance.contexte.nR
    nRc = instance.contexte.nRc
    R = instance.contexte.R
    E=instance.contexte.E
    Cb = instance.contexte.Cb
    Cu = instance.contexte.Cu
    W = instance.contexte.W
    P = instance.contexte.P
    w = instance.contexte.w
    p =instance.contexte.p

    G,DG,Vertices, Edges, Arcs=graph_conctruction(instance, Ebis)

    #Ensemble des unités à couvrir.
    Ub=calcul_unites_couverte_bases(instance)
    U=[(u,t) for u in 1:instance.contexte.nU, t in 1:instance.contexte.nT if (u,t) ∉  Ub]


    # Variables
    if relaxed
        @variables(modele, begin
            x[1:nH, 1:nE]>= 0
            y[1:nR, 1:nE]>= 0
            z[U]>= 0
            l[Arcs, 1:nC]>= 0
            f[Arcs, Vertices[2:size(Vertices,1)]] >= 0
            end)
    else
        @variables(modele, begin
            x[1:nH, 1:nE], Bin
            y[1:nR, 1:nE], Bin
            z[U] >= 0
            f[Arcs, Vertices[2:size(Vertices,1)]] >= 0
            l[Arcs, 1:nC], Bin
            end)
    end

        pos_base = instance.contexte.Eb[1]

        if relaxed
            @constraint(modele, c_r1[h=1:nH, pos=1:nE], x[h,pos]<=1)
            @constraint(modele, c_r2[r=1:nR, pos=1:nE], y[r,pos]<=1)
            @constraint(modele, c_r3[u in U], z[u]<=1)
            @constraint(modele, c_r4[a in Arcs,v in Vertices[2:size(Vertices,1)]], f[a,v]<=1)
            @constraint(modele, c_r5[a in Arcs, c=1:nC], l[a,c]<=1)
        end

    # Contraintes
    @constraint(modele, c_01[h=1:nH, pos in Ebisbar], x[h,pos]==0)

    @constraint(modele, c_02[r=1:nR, pos in Ebisbar], y[r,pos] == 0)

    @constraint(modele, c_03[h=1:nH, pos in Ebis], x[h,pos]  <= sum(y[r,pos] for r in 1:nR))

    @constraint(modele, c_1[h=1:nH], sum(x[h,pos] for pos in Ebis) <= 1)

    @constraint(modele, c_2[pos in Ebis], sum(x[h,pos] for h in 1:nH) <= 1)

    @constraint(modele, c_3[r=1:nR], sum(y[r,pos] for pos in Ebis) <= 1)

    @constraint(modele, c_4[pos in Ebis], sum(w[r] * y[r,pos] for r in 1:nR) <= sum(W[h] * x[h,pos] for h in 1:nH))

    @constraint(modele, c_5[pos in Ebis], sum(p[r] * y[r,pos] for r in 1:nR) <= sum(P[h] * x[h,pos] for h in 1:nH))

    @constraint(modele, c_6[u in U], z[u] <= sum(y[r,pos] for c in 1:length(Cu[u[1]]), r in R[Cu[u[1]][c]], pos in E1[u[1]][c][r][u[2]]) )

    @constraint(modele, c_r33[u in U], z[u]<=1)
    # @constraint(modele, c_7[b=1:nB], sum(y[r,pos] for c in Cb[b], r in R[c], pos in E2[b][r]) >= 1)
#

    # #Flow conservation f
    #
    for v in Vertices[2:size(Vertices,1)]
            for i in 1:nv(DG)
                if i!=v && i!=1
                    @constraint(modele, sum(f[(i,j),v] for j in outneighbors(DG, i)) -
                    sum(f[(j,i),v]  for j in inneighbors(DG, i))  ==0)
                end
                if i==v
                    @constraint(modele, sum(f[(i,j),v]  for j in outneighbors(DG, i)) -
                    sum( f[(j,i),v] for j in inneighbors(DG, i))  == - sum(x[h,Ebis[v-1]] for h in 1:nH))
                end
                if i==1
                    @constraint(modele, sum(f[(i,j),v]  for j in outneighbors(DG, i)) -
                    sum(f[(j,i),v] for j in inneighbors(DG, i))  == sum(x[h,Ebis[v-1]] for h in 1:nH))
                end
            end
        end


# # Defining l variables
for a in Arcs
        for c in 1:nC
            #finding the relays of type c that can reach j from i
            R_cj=Int64[]
            if a[1]!=1 && a[2]!=1
                position = Ebis[a[1]-1]
                for r in R[c]
                    if Ebis[a[2]-1] in E3[position][r]
                        push!(R_cj, r)
                    end
                end
                @constraint(modele, l[a,c] <=
                sum(y[r,Ebis[a[1]-1]] for r in R_cj))

                @constraint(modele, l[a,c] <=
                sum(y[r,Ebis[a[2]-1]] for r in R[c]))

            end

            if a[2]==1
                if c in instance.contexte.Cb[1]
                    @constraint(modele, l[a,c] == sum(y[r,Ebis[a[1]-1]] for r in R[c]
                    if Ebis[a[1]-1] in E2[1][r]))
                else
                    @constraint(modele, l[a,c] == 0)
                end
            end

            if a[1]==1
                if c in instance.contexte.Cb[1]
                    if  Ebis[a[2]-1] in E4[c][1]
                        @constraint(modele, l[a,c] == sum(y[r,Ebis[a[2]-1]] for r in R[c]))
                    else
                        @constraint(modele, l[a,c] == 0)
                    end
                else
                    @constraint(modele, l[a,c] == 0)
                end
            end

        end
end

#Symétrie forte
@constraint(modele,c_sym[a in Arcs, c=1:nC], l[(a[2],a[1]),c] == l[(a[1],a[2]), c])

#Flow interdiction
for a in Arcs
    for v in Vertices[2:size(Vertices,1)]
        @constraint(modele, f[a,v] <= sum(l[a,c] for c in 1:nC))
    end
end

for a in Arcs
    for v in Vertices[2:size(Vertices,1)]
        @constraint(modele, f[a,v] <= sum(x[h,Ebis[v-1]] for h in 1:nH))
    end
end

for v in Vertices[2:size(Vertices,1)]
    @constraint(modele, sum(f[(i,1),v]  for i in inneighbors(DG, 1)) == 0)
    @constraint(modele, sum(f[(v,i),v]  for i in outneighbors(DG, v)) == 0)
end

    if ajout_R1 == true
        @constraint(modele, c_R1[pos in Ebis, c=1:nC], sum(y[r,pos] for r in R[c]) <= sum(x[h, pos] for h in 1:nH))
    end

    if ajout_R2 == true
        @constraint(modele, c_R22[pos in Ebis, r=1:nR], y[r, pos] <= sum(x[h, pos] for h in 1:nH))
    end

    #Relays cover inequalities
    for h in 1:nH
        Rr=[[r for r in R[c]] for c in 1:nC]
        for c in 1:nC
            sort!(Rr[c], by = r -> w[r])
        end
            weight=0
            Rcover=Int[]
            cpt=1
            while weight <= W[h] && cpt<=nC
                weight+=w[Rr[cpt][1]]
                push!(Rcover, Rr[cpt][1])
                deleteat!(Rr[cpt], 1)
                # println(weight)
                cpt+=1
            end

            if  weight > W[h] && Rcover!=[]
                capacity=W[h]
                # println("Cover : $Rcover")
                # println("weight : $weight", "weight capacitty: $capacity")
                for pos in Ebis
                    @constraint(modele, sum(y[r, pos] for r in Rcover) <= length(Rcover)- x[h, pos])
                end
            end
        end


        for h in 1:nH
            Rr=[[r for r in R[c]] for c in 1:nC]
            for c in 1:nC
                sort!(Rr[c], by = r ->p[r])
            end
                weight=0
                Rcover=Int[]
                cpt=1
                while weight <= P[h] && cpt<=nC
                    weight+=p[Rr[cpt][1]]
                    push!(Rcover, Rr[cpt][1])
                    deleteat!(Rr[cpt], 1)
                    # println(weight)
                    cpt+=1
                end

                if  weight > P[h] && Rcover!=[]
                    capacity=P[h]
                    # println("Cover : $Rcover")
                    # println("weight : $weight ", "power capacitty: $capacity")
                    for pos in Ebis
                        @constraint(modele, sum(y[r, pos] for r in Rcover) <= length(Rcover)- x[h, pos])
                    end
                end
            end

    # end


    epsilone=0.000001
    @objective(modele, Max, sum(z[u] for u in U))


    # lazy_callback_called = false
    # function my_lazy_callback_function(cb_data)
    #     lazy_callback_called = true
    #     x_vals = callback_value.(Ref(cb_data), x)
    #     y_vals = callback_value.(Ref(cb_data), y)
    #
    #     # #Determining a relays cover thanks to
    #     for v in Vertices[2:size(Vertices,1)]
    #         weight=0
    #         Rcover=Int[]
    #         Rr=[[r for r in R[c]] for c in 1:nC]
    #
    #         for c in 1:nC
    #             sort!(Rr[c], by = x -> instance.contexte.s[x]/p[x])
    #         end
    #
    #         if sum(x_vals[h,Ebis[v-1]] for h in 1:nH) > 0
    #             cpt=1
    #             while weight <= sum(x_vals[h,Ebis[v-1]]*P[h] for h in 1:nH) && cpt <= nC
    #                 # println("hello ")
    #                 weight+=p[Rr[cpt][1]]
    #                 push!(Rcover, Rr[cpt][1])
    #                 deleteat!(Rr[cpt], 1)
    #                 # println(weight)
    #                 cpt+=1
    #             end
    #
    #             # adding the cut constraint
    #             capacity=sum(x_vals[h,Ebis[v-1]]*P[h] for h in 1:nH)
    #             # println("power : $weight, capacity: $capacity")
    #             println("Cover : $Rcover")
    #
    #             if weight > capacity && Rcover!=[]
    #                 cut_val=sum(y_vals[r,Ebis[v-1]] for r in Rcover)
    #                 Rhs=length(Rcover)-1
    #                 println("cut power value: $cut_val, Rhs: $Rhs")
    #                 if cut_val > Rhs
    #                         con = @build_constraint(
    #                         sum(y[r,Ebis[v-1]] for r in Rcover) <=  length(Rcover)-1
    #                         )
    #                         println("Adding $(con)")
    #                         #MOI.submit(modele, MOI.UserCut(cb_data), con)
    #                         MOI.submit(modele, MOI.LazyConstraint(cb_data), con)
    #             end
    #         end
    #     end
    # end
    #
    #     # #Determining a relays cover
    #     for v in Vertices[2:size(Vertices,1)]
    #         weight=0
    #         Rcover=Int[]
    #         Rr=[[r for r in R[c]] for c in 1:nC]
    #
    #         for c in 1:nC
    #             sort!(Rr[c], by = x -> instance.contexte.s[x]/w[x])
    #         end
    #
    #         if sum(x_vals[h,Ebis[v-1]] for h in 1:nH) > 0
    #             cpt=1
    #             while weight <= sum(x_vals[h,Ebis[v-1]]*W[h] for h in 1:nH) && cpt <= nC
    #                 # println("hello ")
    #                 weight+=w[Rr[cpt][1]]
    #                 push!(Rcover, Rr[cpt][1])
    #                 deleteat!(Rr[cpt], 1)
    #                 # println(weight)
    #                 cpt+=1
    #             end
    #
    #             # adding the cut constraint
    #             capacity=sum(x_vals[h,Ebis[v-1]]*W[h] for h in 1:nH)
    #             # println("weight : $weight, capacity: $capacity")
    #             println("Cover : $Rcover")
    #
    #             if weight > capacity && Rcover!=[]
    #                 cut_val=sum(y_vals[r,Ebis[v-1]] for r in Rcover)
    #                 Rhs=length(Rcover)-1
    #                 println("cut weight value: $cut_val, Rhs: $Rhs")
    #                 if cut_val > Rhs
    #                         con = @build_constraint(
    #                         sum(y[r,Ebis[v-1]] for r in Rcover) <=  length(Rcover)-1
    #                         )
    #                         println("Adding $(con)")
    #                         #MOI.submit(modele, MOI.UserCut(cb_data), con)
    #                         MOI.submit(modele, MOI.LazyConstraint(cb_data), con)
    #                 end
    #             end
    #         end
    #     end
    # end

###CUt callback
        callback_called = false
        function my_callback_function(cb_data)
            callback_called = true
        x_vals = callback_value.(Ref(cb_data), x)
        y_vals = callback_value.(Ref(cb_data), y)
#Power cover inequalities
        for v in Vertices[2:size(Vertices,1)]
            for h in 1:nH
                if x_vals[h,Ebis[v-1]] >0
                    weight=0
                    Rcover=Int[]
                    Rr=[[r for r in R[c]] for c in 1:nC]
                    for c in 1:nC
                        sort!(Rr[c], by = r -> y_vals[r,Ebis[v-1]]*instance.contexte.s[r]/p[r], rev=true)
                    end
                    cpt=1
                    while weight <= P[h] && cpt <= nC
                        weight+=p[Rr[cpt][1]]
                        push!(Rcover, Rr[cpt][1])
                        deleteat!(Rr[cpt], 1)
                        # println(weight)
                        cpt+=1
                    end
                    # println("Cover : $Rcover")
                    if  weight > P[h] && Rcover!=[]
                        cut_val=sum(y_vals[r,Ebis[v-1]] for r in Rcover)
                        Rhs=length(Rcover)-x_vals[h,Ebis[v-1]]
                        # println("cut power value: $cut_val, Rhs: $Rhs")
                        if cut_val > Rhs
                                con = @build_constraint(
                                sum(y[r,Ebis[v-1]] for r in Rcover) <=  length(Rcover)-x[h,Ebis[v-1]]
                                )
                                # println("Adding $(con)")
                                #MOI.submit(modele, MOI.UserCut(cb_data), con)
                                MOI.submit(modele, MOI.UserCut(cb_data), con)
                        end
                    end
                end
            end
        end

#weight cover inequalities
        for v in Vertices[2:size(Vertices,1)]
            for h in 1:nH
                if x_vals[h,Ebis[v-1]]>0
                    weight=0
                    Rcover=Int[]
                    Rr=[[r for r in R[c]] for c in 1:nC]
                    for c in 1:nC
                        sort!(Rr[c], by = r -> y_vals[r,Ebis[v-1]]*instance.contexte.s[r]/p[r], rev=true)
                    end
                    cpt=1
                    while weight <= W[h] && cpt <= nC
                        weight+=w[Rr[cpt][1]]
                        push!(Rcover, Rr[cpt][1])
                        deleteat!(Rr[cpt], 1)
                        # println(weight)
                        cpt+=1
                    end
                    # println("Cover : $Rcover")
                    if  weight > W[h] && Rcover!=[]
                        cut_val=sum(y_vals[r,Ebis[v-1]] for r in Rcover)
                        Rhs=length(Rcover)-x_vals[h,Ebis[v-1]]
                        # println("cut power value: $cut_val, Rhs: $Rhs")
                        if cut_val > Rhs
                                con = @build_constraint(
                                sum(y[r,Ebis[v-1]] for r in Rcover) <=  length(Rcover)-x[h,Ebis[v-1]]
                                )
                                println("Adding $(con)")
                                #MOI.submit(modele, MOI.UserCut(cb_data), con)
                                MOI.submit(modele, MOI.UserCut(cb_data), con)
                        end
                    end
                end
            end
        end

    end

    if relaxed==false
        MOI.set(modele, MOI.UserCutCallback(), my_callback_function)
        # MOI.set(modele, MOI.LazyConstraintCallback(), my_lazy_callback_function)

    end

    return modele
end

function lecture_resultats_modele(instance, E1, E2, modele, Vertices, Arcs, Edges, Ebis, relaxed)
    Ub=calcul_unites_couverte_bases(instance)
    U=[(u,t) for u in 1:instance.contexte.nU, t in 1:instance.contexte.nT if (u,t) ∉  Ub]
    if relaxed
        x = value.(modele[:x])
        z = value.(modele[:z])
        y = value.(modele[:y])
        l = value.(modele[:l])

        for h in 1:instance.contexte.nH
            for pos in 1:instance.contexte.nE
                xx=x[h,pos]
                if xx !=0
                    println("pos, h: $pos,$h, x: $xx")
                end
            end
        end

        for r in 1:instance.contexte.nR
            for pos in 1:instance.contexte.nE
                yy=y[r,pos]
                if yy !=0
                    println("pos, r: $pos,$r, y: $yy")
                end
            end
        end

    for a in Arcs
        for c in 1:instance.contexte.nC
            s=l[a,c]
                if s !=0
                    println("a,c: $a,$c, l: $s")
                end
        end
    end

    for u in U
        if z[u]!=0
            zz=z[u]
            println("u,t: $u, z: $zz")
    end
end
        return true
    else
        epsilone=0.001
        E=instance.contexte.E
        nC=instance.contexte.nC
        pos_base = instance.contexte.Eb[1]

        x = round.(Int, value.(modele[:x]))
        y = round.(Int, value.(modele[:y]))
        f = round.(Int, value.(modele[:f]))
        l = round.(Int, value.(modele[:l]))

        println("\n========== Affichage des résultats ==========")
        z1 = trunc(Int, objective_value(modele))
        println("\n --> Nombre total d'unités couvertes par la base: ", size(Ub,1))
        println("\n --> Nombre total d'unités couvertes par les HAPS/nombre d'unité non couvertes par les bases: ", z1, "/", size(U,1), " (",trunc((z1/(size(U,1)))*100, digits=2),"%)")
        println("\n --> Nombre total d'unités couvertes/ nombre total d'unité: ", size(Ub,1)+z1, "/", instance.contexte.nU*instance.contexte.nT, " (",trunc(((size(Ub,1)+z1)/(instance.contexte.nU*instance.contexte.nT))*100, digits=2),"%)")

        deploiement_haps = [-1 for i in 1:instance.contexte.nH] # deploiement_haps[h] : position sur laquelle est déployé le HAPS h (si déployé)
        placement_relais = [-1 for i in 1:instance.contexte.nR] # placement_relais[r] : HAPS sur lequel est placé le relais r (si placé)
        for h in 1:instance.contexte.nH
            deploye = false
            for pos in 1:instance.contexte.nE
                if x[h,pos] == 1
                    position = instance.contexte.E[pos]
                    println("\n - HAPS $h déployé à la position ($(position.x), $(position.y))")
                    deploiement_haps[h] = pos
                    for c in 1:instance.contexte.nC
                        for r in instance.contexte.R[c]
                            if y[r,pos] == 1
                                println("\n     - Relais $r (type $c) placé dans ce HAPS")
                                placement_relais[r] = h
                            end
                        end
                    end
                end
            end
        end

        z = round.(Int, value.(modele[:z])) # On arrondit car Cbc renvoie parfois 0.999999... (tolérance)
        unites_couvertes = Vector{Vector{Vector{Vector{Int}}}}(undef, instance.contexte.nU) # unites_couvertes[u][t][c] : ensemble des HAPS couvrant l'unité u au temps t pour le type de communication c (dont dispose l'unité u)
        nb_unites_non_couvertes = [0 for i in 1:instance.contexte.nT]
        for u in 1:instance.contexte.nU
            unites_couvertes[u] = Vector{Vector{Vector{Int}}}(undef, instance.contexte.nT)
            for t in 1:instance.contexte.nT
                unites_couvertes[u][t] = [Int[] for c in 1:instance.contexte.nC]
                if (u,t) in U
                    if z[(u,t)] == 0 # L'unité u est couverte au temps t
                        for c in 1:length(instance.contexte.Cu[u])
                            # On cherche par quel(s) types de communication
                            for h in 1:instance.contexte.nH
                                # Par quel HAPS
                                for r in instance.contexte.R[instance.contexte.Cu[u][c]]
                                    # Quel relais
                                    for pos in E1[u][c][r][t]
                                        # Positions potentiellement couvrantes si relais r placé
                                        if x[h,pos] == 1 && y[r,pos] == 1
                                            # Si HAPS déployé ET relais déployé alors unité couverte par ce HAPS et pour ce type de communication
                                            push!(unites_couvertes[u][t][instance.contexte.Cu[u][c]], h)
                                            break
                                        end
                                    end
                                end
                            end
                        end
                elseif z[(u,t)] == 1
                    nb_unites_non_couvertes[t] += 1
                end
            end
            end
        end

        bases_couvertes = Vector{Vector{Vector{Int}}}(undef, instance.contexte.nB) # unites_couvertes[b][c] : ensemble des HAPS couvrant l'unité b pour le type de communication c (dont dispose l'unité u)
        for b in 1:instance.contexte.nB
            bases_couvertes[b] = [Int[] for c in 1:instance.contexte.nC]
            for c in instance.contexte.Cb[b]
                for h in 1:instance.contexte.nH
                # Par quel HAPS
                    for r in instance.contexte.R[c]
                        # Quel relais
                        for pos in E2[b][r]
                            # Positions potentiellement couvrantes si relais r placé
                            if x[h,pos] == 1 && y[r,pos] == 1
                                # Si HAPS déployé ET relais déployé alors base couverte par ce HAPS et pour ce type de communication
                                push!(bases_couvertes[b][c], h)
                                break
                            end
                        end
                   end
                end
            end
        end

        println("\n========== Affichage des résultats ==========")
        #flow variables
            for a in Arcs
                for v in  Vertices[2:size(Vertices,1)]
                        flow=f[a,v]
                        if flow !=0
                            if a[1]==1
                                e1=pos_base
                                e2=E[Ebis[a[2]-1]]

                                h1=10000
                                for h in 1:instance.contexte.nH
                                        if x[h,Ebis[a[2]-1]] == 1
                                            h2=h
                                        end
                                end

                            end
                            if a[2]==1
                                e2=pos_base
                                e1=E[Ebis[a[1]-1]]

                                h2=10000
                                for h in 1:instance.contexte.nH
                                        if x[h,Ebis[a[1]-1]] == 1
                                            h1=h
                                        end
                                end


                            end
                            if a[1]!=1 && a[2]!=1
                                e2=E[Ebis[a[2]-1]]
                                e1=E[Ebis[a[1]-1]]

                                for h in 1:instance.contexte.nH
                                    for h_bis in 1:instance.contexte.nH
                                        if x[h,Ebis[a[1]-1]] == 1 && x[h_bis, Ebis[a[2]-1]] == 1
                                            h1=h
                                            h2=h_bis
                                        end
                                    end
                                end

                            end
                        println("\n - flow f $v arc $a, ($h1, $h2), pos ($e1, $e2) : $flow")
                        end
                end
        end


    transmission_haps = Vector{Vector{Vector{Int}}}(undef, instance.contexte.nH) # transmission_haps[h][c] : HAPS vers lequel le HAPS h peut transmettre via le type de communication c
    for h in 1:instance.contexte.nH
        transmission_haps[h] = [Int[] for c in 1:instance.contexte.nC]
    end

    R=instance.contexte.R
    for a in Arcs
        for c in 1:instance.contexte.nC
            if a[1]!=1 && a[2]!=1
                pos=Ebis[a[1]-1]
                pos_bis=Ebis[a[2]-1]
                if sum(f[a,v] for v in Vertices[2:size(Vertices,1)])!=0 &&
                    sum(y[r,pos] for r in R[c]) !=0 &&
                    sum(y[r,pos_bis] for r in R[c]) !=0
                        for h in 1:instance.contexte.nH
                            for h_bis in 1:instance.contexte.nH
                                if x[h,pos] == 1 && x[h_bis, pos_bis] == 1
                                    push!(transmission_haps[h][c], h_bis)
                                    println((h,h_bis,c))
                                end
                            end
                        end
                    end
                end
            end
    end


        solution = Solution(1, z1, nb_unites_non_couvertes, deploiement_haps, placement_relais, unites_couvertes, bases_couvertes, transmission_haps)
    end
    return solution
end

function positions_non_couvrantes(instance)
    seuil_max_relais = maximum(instance.contexte.s)
    E = instance.contexte.E
    total = 0
    for position in 1:instance.contexte.nE
        position_active = false
        t = 1
        while t <= instance.contexte.nT && position_active == false
             u = 1
             while u <= instance.contexte.nU && position_active == false
                pos_unite = instance.scenario.deplacements[u][t]
                if distance(instance, pos_unite, E[position]) <= seuil_max_relais
                    position_active = true # Position à portée d'au moins une unité
                end
                u += 1
             end
            t += 1
        end
        b = 1
        while b <= instance.contexte.nB && position_active == false
            pos_base = instance.contexte.Eb[b]
            if distance(instance, pos_base, E[position]) <= seuil_max_relais
                position_active = true # Position à portée d'au moins une base
            end
            b += 1
        end
        if position_active == false
            total += 1
            instance.contexte.Ebool[position] = false
        end
    end
    total_pourcentage = trunc((total/instance.contexte.nE)*100, digits=2)
    println("\n - Nombre de positions désactivées : $total/$(instance.contexte.nE) ($(total_pourcentage)%)")
end



function distance(instance, p1, p2)
    A = instance.contexte.A
    return sqrt((p2.x - p1.x)^2 + (p2.y - p1.y)^2 + A^2)
end

function calcul_positions_unites(instance)
    E1 = Vector{Vector{Vector{Vector{Vector{Int}}}}}(undef, instance.contexte.nU)
    # E1[u][r][t] : ensemble des positions à portée de l'unité u de type c pour le relais r
    for u in 1:instance.contexte.nU
        E1[u] = Vector{Vector{Vector{Vector{Int}}}}(undef,length(instance.contexte.Cu[u]))
            for c in 1:length(instance.contexte.Cu[u])
                E1[u][c] = Vector{Vector{Vector{Int}}}(undef, instance.contexte.nR)
                for r in 1:instance.contexte.nR
                    E1[u][c][r] = Vector{Vector{Int}}(undef, instance.contexte.nT)
                    for t in 1:instance.contexte.nT
                        E1[u][c][r][t] = Int[]
                        for index_position in 1:instance.contexte.nE
                            position = instance.contexte.E[index_position]
                            pos_unite = instance.scenario.deplacements[u][t]
                            if distance(instance, pos_unite, position) <= min(instance.contexte.s[r],instance.contexte.Pu[u][c])
                                if !(index_position in E1[u][c][r][t])
                                    push!(E1[u][c][r][t], index_position)
                                end
                            end
                        end
                    end
                end
            end
    end
    return E1
end

function calcul_positions_bases(instance)
    E2 = Vector{Vector{Vector{Int}}}(undef, instance.contexte.nB)
    # E2[b][r] : ensemble des positions à portée de la base b pour le relais r
    for b in 1:instance.contexte.nB
        E2[b] = Vector{Vector{Int}}(undef, instance.contexte.nR)
        for r in 1:instance.contexte.nR
            E2[b][r] = Int[]
            for index_position in 1:instance.contexte.nE
                position = instance.contexte.E[index_position]
                pos_base = instance.contexte.Eb[b]
                if distance(instance, pos_base, position) <= instance.contexte.s[r]
                    if !(index_position in E2[b][r])
                        push!(E2[b][r], index_position)
                    end
                end
            end
        end
    end
    return E2
end

function graph_conctruction(instance, Ebis)
    E=instance.contexte.E
    Vertices=Int64[]
    Edges=Tuple{Int64,Int64}[]
    Arcs=Tuple{Int64,Int64}[]

    G = complete_graph(size(Ebis,1)+1)

    push!(Vertices, 1)
    cmp=2
    for e in Ebis
        push!(Vertices, cmp)
        cmp+=1
    end


    for e in edges(G)
        u, v = src(e), dst(e)
        push!(Edges, (u,v))
    end

    #remove the uninteressant edges
    pos_base = instance.contexte.Eb[1]
    seuil_max_relais = maximum(instance.contexte.s)

    for a in edges(G)
            u, v = src(a), dst(a)
            if u==1
                e1=pos_base
                e2=E[Ebis[v-1]]
            end
            if v==1
                e1=E[Ebis[u-1]]
                e2=pos_base
            end
            if u!=1 && v!=1
                e1=E[Ebis[u-1]]
                e2=E[Ebis[v-1]]
            end

            if distance(instance, e1, e2) > seuil_max_relais
                rem_edge!(G, a)
            end
    end

    DG=DiGraph(G)

    for a in edges(DG)
            u, v = src(a), dst(a)
            push!(Arcs, (u,v))
    end

Edges=Tuple{Int64,Int64}[]

for e in edges(G)
    u, v = src(e), dst(e)
    push!(Edges, (u,v))
end

return G,DG,Vertices, Edges, Arcs
end

function edge_(a,Edges)
    for e in Edges
        if (a[1]==e[1] && a[2]==e[2]) ||(a[1]==e[2] && a[2]==e[1])
            return e
        end
    end
end


function distance_air_air(instance, p1, p2)
    return sqrt((p2.x - p1.x)^2 + (p2.y - p1.y)^2)
end


function calcul_positions_relais(instance, Ebis)
    E3 = Vector{Vector{Vector{Int}}}(undef, instance.contexte.nE)
    # E3[e][r] : ensemble des positions à portée de la position e pour le relais r
    for position in Ebis
        E3[position] = Vector{Vector{Int}}(undef, instance.contexte.nR)
        for r in 1:instance.contexte.nR
            E3[position][r] = Int[]
            for position_bis in 1:instance.contexte.nE
                p1 = instance.contexte.E[position]
                p2 = instance.contexte.E[position_bis]
                position_active = instance.contexte.Ebool[position]
                position_bis_active = instance.contexte.Ebool[position_bis]
                if distance_air_air(instance, p1, p2) <= instance.contexte.s[r] && position != position_bis && position_active == true && position_bis_active == true
                    push!(E3[position][r], position_bis)
                end
            end
        end
    end
    return E3
end

function calcul_positions_porte_bases(instance)
    E4 =   Vector{Vector{Vector{Int}}}(undef, instance.contexte.nC)
    for c in 1:instance.contexte.nC
        seuil_max_relais = maximum([instance.contexte.s[r] for r in instance.contexte.R[c]])
        E4[c] =   Vector{Vector{Int}}(undef, instance.contexte.nB)
        for b in 1:instance.contexte.nB
            E4[c][b] = Int[]
                for index_position in 1:instance.contexte.nE
                    position = instance.contexte.E[index_position]
                    pos_base = instance.contexte.Eb[b]
                    if distance(instance, pos_base, position) <= seuil_max_relais
                            push!(E4[c][b], index_position)
                    end
                end
        end
    end
        #println("E4: $E4")
return E4
end

function calcul_unites_couverte_bases(instance)
    Ub=Tuple{Int64,Int64}[]
    sr=Int[]
    for b in 1:instance.contexte.nB
        pos_base = instance.contexte.Eb[b]
        for t in 1:instance.contexte.nT
             for u in 1:instance.contexte.nU
                 for c in 1:length(instance.contexte.Cu[u])
                    seuil_max_relais = maximum(instance.contexte.s[r] for r in instance.contexte.R[instance.contexte.Cu[c][1]])
                    pos_unite = instance.scenario.deplacements[u][t]
                    if distance(instance, pos_unite,pos_base) <= min(seuil_max_relais, instance.contexte.Pu[u][c])
                        push!(Ub,(u,t))
                        # d=distance(instance, pos_unite,pos_base)
                        # println("u: $u, t:  $t, distance to base: $d")
                    end
                end
            end
        end
    end
    # println(Ub)
    return Ub
end
