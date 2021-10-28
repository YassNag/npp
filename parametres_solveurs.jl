function gestion_parametres_GLPK(modele)
    choix_verbosite = choix_binaire("\n --> Souhaitez-vous obtenir les détails de la résolution (o/n) ? ")
    if choix_verbosite == "o"
        set_optimizer_attribute(modele, "msg_lev", 3)
    end
    return choix_verbosite
end

function gestion_parametres_CBC(modele)
    choix_verbosite = choix_binaire("\n --> Souhaitez-vous obtenir les détails de la résolution (o/n) ? ")
    if choix_verbosite == "n"
        # CBC verbeux par défaut
        set_optimizer_attribute(modele, "logLevel", 0)
    end
    return choix_verbosite
end

function gestion_parametres_SCIP(modele)
    choix_verbosite = choix_binaire("\n --> Souhaitez-vous obtenir les détails de la résolution (o/n) ? ")
    if choix_verbosite == "n"
        # SCIP verbeux par défaut
        set_optimizer_attribute(modele, "display/verblevel", 0)
    end
    return choix_verbosite
end

function gestion_parametres_CPLEX(modele)
    choix_verbosite = choix_binaire("\n --> Souhaitez-vous obtenir les détails de la résolution (o/n) ? ")
    if choix_verbosite == "n"
        # CPLEX verbeux par défaut
        set_optimizer_attribute(modele, "CPX_PARAM_SCRIND", 0)
    end
    return choix_verbosite
end

function gestion_parametres_Gurobi(modele)
    choix_verbosite = choix_binaire("\n --> Souhaitez-vous obtenir les détails de la résolution (o/n) ? ")
    if choix_verbosite == "n"
        # Gurobi verbeux par défaut
        set_optimizer_attribute(modele, "OutputFlag", 0)    
    end
    return choix_verbosite
end

function gestion_parametres_Mosek(modele)
    choix_verbosite = choix_binaire("\n --> Souhaitez-vous obtenir les détails de la résolution (o/n) ? ")
    if choix_verbosite == "n"
        # Mosek verbeux par défaut
        set_optimizer_attribute(modele, "LOG", 0)    
    end
    return choix_verbosite
end
