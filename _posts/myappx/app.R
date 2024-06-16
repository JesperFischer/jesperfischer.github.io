
pacman::p_load(
  ordbetareg,
  dplyr,
  ggplot2,
  haven,
  shinyWidgets,
  shiny,
  tidyr,
  stringr,
  extraDistr,
  ggh4x
)
# 
# input = data.frame(beta = 1, lapse = 0.05, alpha = 5, conf_int = 1, conf_beta = -5, conf_prec = 10,conf_shift = -2,
#                    rt_int = 2, rt_beta = 2, rt_sd = 1, rt_ndt = 1,rt_shift = 2,conf_cutzero = -1, conf_cutone = 1, id = 1)


# 
# qq = raw_hrd %>% filter(participant_id %in% unique(raw_hrd$participant_id)[170] & Modality == "Extero", session == 1) %>% 
#   mutate(DecisionRT = as.numeric(DecisionRT), Confidence = as.numeric(Confidence), Decision = ifelse(Decision == "Less",0,ifelse(Decision == "More",1,NA))) %>%
#   select(Alpha,DecisionRT,Confidence,Decision)
# 
# qq %>% pivot_longer(cols = c(DecisionRT, Confidence, Decision)) %>% ggplot(aes(x = Alpha, y = value))+
#   geom_point()+geom_smooth()+facet_wrap(~name, ncol = 1, scales = "free")

shiny::shinyApp(
  # Define UI
  ui <- navbarPage(
    "Parameters",
    tabPanel("JPDM", 
             sidebarLayout(
               sidebarPanel(
                 column(width = 4,
                        sliderInput("beta2", "slope", min = -3, max = 3, value = 0.1, step = 0.1),
                        sliderInput("lapse2", "lapse", min = -5, max = 0, value = -4, step = 0.01),
                        sliderInput("alpha2", "alpha", min = -40, max = 40, value = 0, step = 0.1)
                 ),
                 column(width = 4,
                        sliderInput("conf_int2", "conf_int", min = -3, max = 3, value = -1, step = 0.1),
                        sliderInput("conf_beta2", "conf_beta", min = -20, max = 5, value = -5, step = 0.1),
                        sliderInput("conf_shift2", "conf_shift", min = -10, max = 10, value = 0, step = 0.1),
                 ),
                 column(width = 4,
                        sliderInput("rt_int2", "rt_int", min = -3, max = 3, value = -0.5, step = 0.01),
                        sliderInput("rt_beta2", "rt_beta", min = -1, max = 20, value = 5, step = 0.1),
                        sliderInput("rt_shift2", "rt_shift", min = -10, max = 10, value = 0, step = 0.1),
                        sliderInput("rt_ndt2", "rt_ndt", min = 0, max = 2, value = 0.2, step = 0.1),
                 )
               ),
               mainPanel(
                 plotOutput("Relationships", height = "700px")
               )
             )
    ),
    
    tabPanel("Histogram / densities",
             sidebarLayout(
               sidebarPanel(
                 column(width = 4,
                        sliderInput("beta", "slope", min = -3, max = 3, value = 0.1, step = 0.1),
                        sliderInput("lapse", "lapse", min = -5, max = 0, value = -4, step = 0.01),
                        sliderInput("alpha", "alpha", min = -50, max = 50, value = 0, step = 0.1),
                 ),
                 column(width = 4,
                        sliderInput("conf_int", "conf_int", min = -3, max = 3, value = 0, step = 0.1),
                        sliderInput("conf_beta", "conf_beta", min = -20, max = 5, value = -5, step = 0.1),
                        sliderInput("conf_shift", "conf_shift", min = -10, max = 10, value = 0, step = 0.1),
                        sliderInput("conf_cutzero", "conf_cutzero", min = -3, max = 0, value = -2, step = 0.1),
                        sliderInput("conf_cutone", "conf_cutone", min = 0, max = 3, value = 2, step = 0.1),
                        sliderInput("conf_prec", "conf_prec", min = -10, max = 10, value = 3, step = 0.1)
                 ),
                 column(width = 4,
                        sliderInput("rt_int", "rt_int", min = -3, max = 3, value = -0.5, step = 0.01),
                        sliderInput("rt_beta", "rt_beta", min = -1, max = 20, value = 5, step = 0.1),
                        sliderInput("rt_shift", "rt_shift", min = -10, max = 10, value = 0, step = 0.1),
                        sliderInput("rt_ndt", "rt_ndt", min = 0, max = 2, value = 0.2, step = 0.1),
                        sliderInput("rt_sd", "rt_sd", min = -10, max = 10, value = 0, step = 0.1)
                 )
               ),
               mainPanel(
                 plotOutput("Histograms", height = "700px")
               )
             )
    )
  ),
  
  # Define server
  server <- function(input, output) {
    
    
    generate_exact_data <- function() {
      data <- data.frame(
        Alpha = c(40.5, -40.5, -20.5, 20.5, 0.5, -0.5, -12.5, -19.5, -0.5, -7.5, 
                  4.5, 11.5, 16.5, 4.5, -2.5, 4.5, -7.5, -9.5, -5.5, -0.5, 
                  -1.5, 6.5, 2.5, 2.5, 6.5, 5.5, 3.5, 3.5, 0.5, 1.5, 
                  9.5, -9.5, -9.5, 8.5, -11.5, 10.5, -12.5, 11.5, -13.5, 12.5, 
                  -14.5, 12.5, -14.5, -14.5, 13.5, -15.5, 13.5, -15.5, 13.5, -15.5, 
                  13.5, -15.5, 13.5, -15.5, 13.5, -15.5, 13.5, -14.5, -16.5, 13.5, 
                  -16.5, 12.5, -16.5, 12.5, -15.5, 12.5, -15.5, -15.5, 12.5, -15.5, 
                  11.5, -14.5, 11.5, -14.5, 11.5, -14.5, 11.5, -13.5, 10.5, -13.5),
        DecisionRT = c(1.4498747, 2.1154120, 3.1157116, 2.4978739, 6.5660700, 3.1482725, 
                       2.2142246, 1.9667107, 5.4664995, 2.4156643, 5.3823501, 1.6156799, 
                       1.5490398, 5.5147322, 1.5145694, 1.7826786, 2.6826405, 2.4816994, 
                       6.2650572, 2.1663445, 3.8827280, 1.5168280, 4.2817310, 2.9984809, 
                       1.8471111, 2.9008777, 3.1318292, 2.7306845, 2.2971363, 2.0477431, 
                       1.7995315, 2.6497006, 2.1661950, 3.6989570, 4.3816510, 2.5979754, 
                       2.3802903, 2.0334436, 1.7664442, 2.6995966, 3.0664553, 2.8493271, 
                       1.9497760, 4.2654037, 1.6478815, 2.0483511, 1.7989167, 1.4467499, 
                       3.0329406, 1.7986011, 1.4154198, 1.7981120, 1.7815026, 1.8146976, 
                       1.3808055, 1.4135229, 1.5648370, 3.8977399, 2.8163894, 1.7157264, 
                       1.8486031, 1.1662610, 2.1822189, 2.1998760, 1.6484192, 1.9315366, 
                       1.7649654, 1.8666350, 1.4829629, 1.6988575, 1.4657493, 1.3996992, 
                       1.1655448, 2.0153748, 1.3982426, 2.7647399, 1.5647009, 1.6481751, 
                       0.9481878, 1.6981365),
        Confidence = c(100, 98, 100, 100, 2, 37, 100, 100, 37, 59, 
                       0, 99, 100, 46, 77, 82, 56, 62, 11, 55, 
                       36, 93, 37, 47, 75, 72, 53, 47, 74, 59, 
                       89, 87, 90, 54, 44, 85, 100, 75, 60, 76, 
                       53, 76, 70, 18, 89, 100, 100, 56, 74, 94, 
                       99, 81, 100, 100, 100, 100, 100, 2, 100, 99, 
                       93, 100, 100, 100, 90, 100, 94, 100, 100, 100, 
                       98, 100, 100, 77, 100, 71, 100, 100, 100, 98) / 100,
        Decision = c(1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 
                     0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 
                     0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 
                     1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 
                     0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 
                     1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 
                     0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 
                     1, 0, 1, 0, 1, 0, 1, 0, 1, 0)
      )
      return(data)
    }
    

    
    
    output$Relationships <- renderPlot({
      
      inv_logit_scaled = function(x){
        return(exp(x)/(1+exp(x)))
      }
      
      jpdm = function(x,beta,lapse,alpha,
                      conf_int,conf_beta,conf_prec,conf_shift,
                      rt_int,rt_beta,rt_sd,rt_ndt,rt_shift, 
                      id = 1){
        
        probab = function(x,lapse,beta,alpha){
          return((inv_logit_scaled(lapse) / 2) + (1 - 2 *  (inv_logit_scaled(lapse) / 2)) * (0.5+0.5 * pracma::erf((x-alpha) / (exp(beta) * sqrt(2)))))
          
        }
        prob = probab(x,lapse,beta,alpha)
        
        conf = inv_logit_scaled(conf_int + conf_beta * (probab(x-conf_shift,lapse,beta,alpha) * (1-probab(x-conf_shift,lapse,beta,alpha))))
        
        rt = exp(rt_int + rt_beta * (probab(x-rt_shift,lapse,beta,alpha) * (1-probab(x-rt_shift,lapse,beta,alpha)))) + rt_ndt
        
        return(data.frame(prob = prob,rt = rt, conf  = conf,x = x, id = id))
      }
      

      data <- jpdm(x = seq(-50,50,0.1), beta = input$beta2, lapse = input$lapse2, alpha = input$alpha2,
                   conf_int = input$conf_int2, conf_beta = input$conf_beta2, conf_prec = input$conf_prec2, input$conf_shift2,
                   rt_int = input$rt_int2, rt_beta = input$rt_beta2, rt_sd = input$rt_sd2, rt_ndt = input$rt_ndt2, input$rt_shift2, id = 1)
      
      
      
      data = data %>% rename(Probability = prob, Reactiontime = rt, Confidence = conf) %>% 
        pivot_longer(cols = c("Probability","Reactiontime","Confidence"), names_to = "decision",values_to = "decision_value")
      
      
      #overlay real data:
      real_data <- generate_exact_data() %>% rename(Probability = Decision, Reactiontime = DecisionRT ) %>% 
        pivot_longer(cols = c("Probability","Reactiontime","Confidence"), names_to = "decision",values_to = "decision_value")
      
      
      plot1 = data %>% 
        ggplot()+
        geom_line(data = data,aes(x = x, y = decision_value))+
        geom_point(data = real_data, aes(x = Alpha, y = decision_value), fill = "lightblue",color = "black", size = 3, shape = 21)+
        facet_wrap(~decision, scales = "free",ncol = 1)+theme_classic()+
        theme(text = element_text(size = 24))
      
      
      plot1 = plot1+
        facetted_pos_scales(
          y = list(
            decision == "Probability" ~ scale_y_continuous(limits = c(0, 1), breaks = c(0,0.2,0.4,0.6,0.8,1.0)),
            decision == "Reactiontime" ~ scale_y_continuous(limits = c(0, 8), breaks = c(0,2,4,6,8)),
            decision == "Confidence" ~ scale_y_continuous(limits = c(0, 1), breaks = c(0,0.2,0.4,0.6,0.8,1.0))
          )
        )
      
      
      plot1
    })
    
    
    output$Histograms <- renderPlot({
      
      inv_logit_scaled = function(x){
        return(exp(x)/(1+exp(x)))
      }
      
      rordbeta <- function(n=100,
                           mu=0.5,
                           phi=1,
                           cutpoints=c(-1,1)) {
        
        if(!all(mu>0 & mu<1) ) {
          
          stop("Please pass a numeric value for mu that is between 0 and 1.")
          
        }
        
        if(!all(phi>0)) {
          
          stop("Please pass a numeric value for phi that is greater than 0.")
          
        }
        
        if(!(length(mu) %in% c(1,n))) {
          
          stop("Please pass a vector for mu that is either length 1 or length N.")
          
        }
        
        if(!(length(phi) %in% c(1,n))) {
          
          stop("Please pass a vector for phi that is either length 1 or length N.")
          
        }
        
        mu_ql <- qlogis(mu)
        
        if(length(mu_ql)==1) {
          mu_ql <- rep(mu_ql, n)
        }
        
        # probabilities for three possible categories (0, proportion, 1)
        low <- 1-plogis(mu_ql - cutpoints[1])
        middle <- plogis(mu_ql - cutpoints[1]) - plogis(mu_ql - cutpoints[2])
        high <- plogis(mu_ql - cutpoints[2])
        
        # we'll assume the same eta was used to generate outcomes
        
        out_beta <- rbeta(n = n,mu * phi, (1 - mu) * phi)
        
        # now determine which one we get for each observation
        outcomes <- sapply(1:n, function(i) {
          
          sample(1:3,size=1,prob=c(low[i],middle[i],high[i]))
          
        })
        
        # now combine binary (0/1) with proportion (beta)
        
        final_out <- sapply(1:n,function(i) {
          if(outcomes[i]==1) {
            return(0)
          } else if(outcomes[i]==2) {
            return(out_beta[i])
          } else {
            return(1)
          }
        })
        
        
        return(final_out)
        
      }
      
      
      jpdm = function(x,beta,lapse,alpha,
                      conf_int,conf_beta,conf_prec,conf_shift,conf_cutzero,conf_cutone,
                      rt_int,rt_beta,rt_sd,rt_ndt,rt_shift, 
                      id = 1){
        
        probab = function(x,lapse,beta,alpha){
          return((inv_logit_scaled(lapse) / 2) + (1 - 2 *  (inv_logit_scaled(lapse) / 2)) * (0.5+0.5 * pracma::erf((x-alpha) / (exp(beta) * sqrt(2)))))
          
        }
        prob = probab(x,lapse,beta,alpha)
        prob = rbinom(length(prob),1,prob)
        
        conf = inv_logit_scaled(conf_int + conf_beta * (probab(x-conf_shift,lapse,beta,alpha) * (1-probab(x-conf_shift,lapse,beta,alpha))))
        
        conf = rordbeta(length(conf),conf,exp(conf_prec),c(conf_cutzero,conf_cutone))
        
        rts = (rt_int + rt_beta * (probab(x-rt_shift,lapse,beta,alpha) * (1-probab(x-rt_shift,lapse,beta,alpha))))
        rt = rlnorm(length(rts),rts, exp(rt_sd)) + rt_ndt
        
        return(data.frame(prob = prob,rt = rt, conf  = conf,x = x, id = id))
      }
      
      
      data <- jpdm(x = seq(-50,50,0.01), beta = input$beta, lapse = input$lapse, alpha = input$alpha,
                   conf_int = input$conf_int, conf_beta = input$conf_beta, conf_prec = input$conf_prec, input$conf_shift,input$conf_cutzero,input$conf_cutone,
                   rt_int = input$rt_int, rt_beta = input$rt_beta, rt_sd = input$rt_sd, rt_ndt = input$rt_ndt, input$rt_shift, id = 1)
      
      
      
      data = data %>% rename(Probability = prob, Reactiontime = rt, Confidence = conf) %>% 
        pivot_longer(cols = c("Probability","Reactiontime","Confidence"), names_to = "decision",values_to = "decision_value")
      
      real_data <- generate_exact_data() %>% rename(Probability = Decision, Reactiontime = DecisionRT ) %>% 
        pivot_longer(cols = c("Probability","Reactiontime","Confidence"), names_to = "decision",values_to = "decision_value")
      
      
      
      plot1 = data %>% 
        ggplot()+
        geom_histogram(data = data,aes(x = decision_value, y = after_stat(density)), position = "identity", col = "black", bins = 30)+
        geom_histogram(data = real_data,aes(x = decision_value, y = after_stat(density)),
                       position = "identity", alpha = 0.5, col = "black", bins = 30, fill = "lightblue")+
        
        facet_wrap(~decision, scales = "free",ncol = 1)+theme_classic()+
        theme(text = element_text(size = 24))
      
      plot1 = plot1+
        facetted_pos_scales(
          x = list(
            decision == "Confidence" ~ scale_x_continuous(limits = c(0, 1), breaks = c(0,0.2,0.4,0.6,0.8,1.0)),
            decision == "Probability" ~ scale_x_continuous(limits = c(0, 1), breaks = c(0,1.0)),
            decision == "Reactiontime" ~ scale_x_continuous(limits = c(0, 10), breaks = c(0,1,2,3,4,5,6,7,8,9,10))
          )
        )
      
      plot1
      
      
      
    })
    
    
    
    
  },
  options = list(height = 700, width = 700)
  
)
