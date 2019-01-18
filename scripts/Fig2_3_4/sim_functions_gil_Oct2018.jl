### Gillespie implementation

### WT
function sim_delay(p::Array{Float64,1},tspan::Tuple{Float64,Float64},u0::Array{Int64,1})
  t::Float64=tspan[1]
  i::Int64=0
  out::Array{Array{Int64,1}} = []
  t_out::Array{Float64,1} = []
  push!(out, u0)
  push!(t_out,t)
  u::Array{Int64,1}=copy(u0)
  rate::Array{Float64,1}=zeros(56)
  while (t<tspan[2])
    i+=1
    # Xist transcription
    rate[1] = p[21]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[7]^p[13]/(u[7]^p[13]+(p[22]*p[14])^p[13]))
    rate[2] = p[21]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[8]^p[13]/(u[8]^p[13]+(p[22]*p[14])^p[13]))
    # cXR production
    rate[3] = p[22]*(1-u[51]^p[5]/(u[51]^p[5]+(p[21]*p[6])^p[5]))
    rate[4] = p[22]*(1-u[49]^p[5]/(u[49]^p[5]+(p[21]*p[6])^p[5]))
    # tXA production
    rate[5] = p[23]*(1-u[52]^p[3]/(u[52]^p[3]+(p[21]*p[4])^p[3]))
    rate[6] = p[23]*(1-u[50]^p[3]/(u[50]^p[3]+(p[21]*p[4])^p[3]))
    # degradation
    rate[7] = 0.1733*u[3]
    rate[8] = 0.1733*u[4]
    rate[9] = 1*u[7]
    rate[10] = 1*u[8]
    rate[11] = 1*u[5]
    rate[12] = 1*u[6]
    rate[13] = 1*u[49]
    rate[14] = 1*u[50]
    rate[15] = 1*u[51]
    rate[16] = 1*u[52]
    ### silencing intermediates on chr1
    rate[17] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[1],0)
    rate[18] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[10],0)
    rate[19] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[11],0)
    rate[20] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[12],0)
    rate[21] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[13],0)
    rate[22] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[14],0)
    rate[23] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[15],0)
    rate[24] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[16],0)
    rate[25] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[17],0)
    rate[26] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[18],0)
    rate[27] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[19],0)
    rate[28] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[20],0)
    rate[29] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[21],0)
    rate[30] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[22],0)
    rate[31] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[23],0)
    rate[32] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[24],0)
    rate[33] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[25],0)
    rate[34] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[26],0)
    rate[35] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[27],0)
    rate[36] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[28],0)
    ### silencing intermediates chr 2
    rate[37] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[2],0)
    rate[38] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[30],0)
    rate[39] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[31],0)
    rate[40] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[32],0)
    rate[41] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[33],0)
    rate[42] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[34],0)
    rate[43] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[35],0)
    rate[44] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[36],0)
    rate[45] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[37],0)
    rate[46] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[38],0)
    rate[47] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[39],0)
    rate[48] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[40],0)
    rate[49] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[41],0)
    rate[50] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[42],0)
    rate[51] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[43],0)
    rate[52] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[44],0)
    rate[53] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[45],0)
    rate[54] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[46],0)
    rate[55] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[47],0)
    rate[56] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[48],0)
	
	# Gillespie
    r = rand(2)
    # determine time point of next reaction
    t_step::Float64 = 1/sum(rate)*log(1/(1-r[1]))
    # decide which reaction will occur
    temp::Float64 = r[2]*sum(rate)
    cum_rate = cumsum(rate)
    # Execute reaction
    if temp<cum_rate[1]
      u[1] += 1
      u[3] += 1
        if p[7]==0
            u[49] += 1
        end
        if p[8]==0
            u[50] += 1
        end
  elseif temp<cum_rate[2]
      u[2] += 1
      u[4] += 1
      if p[7]==0
        u[51] += 1
      end
      if p[8]==0
        u[52] += 1
      end
  elseif temp<cum_rate[3]
        u[8] += 1
  elseif temp<cum_rate[4]
      u[7] += 1
  elseif temp<cum_rate[5]
        u[6] += 1
  elseif temp<cum_rate[6]
        u[5] += 1
  elseif temp<cum_rate[7]
      if maximum([p[7],p[8]])==0
        u[1] -= 1
      end
        u[3] -=1
  elseif temp<cum_rate[8]
      if maximum([p[7],p[8]])==0
        u[2] -= 1
      end
        u[4] -=1
  elseif temp<cum_rate[9]
        u[7] -= 1
  elseif temp<cum_rate[10]
        u[8] -= 1
  elseif temp<cum_rate[11]
        u[5] -= 1
  elseif temp<cum_rate[12]
        u[6] -= 1
  elseif temp<cum_rate[13]
        u[49] -= 1
  elseif temp<cum_rate[14]
        u[50] -= 1
  elseif temp<cum_rate[15]
        u[51] -= 1
  elseif temp<cum_rate[16]
        u[52] -= 1
  elseif temp<cum_rate[17]
      if maximum([p[7],p[8]])>1
        u[10] += 1
        u[1] -= 1
      else
        u[1] -= 1
      end
      p[7]==1 ? u[49] += 1 : 0
      p[8]==1 ? u[50] += 1 : 0
  elseif temp<cum_rate[18]
      if maximum([p[7],p[8]])>2
        u[10] -= 1
        u[11] += 1
      else
        u[10] -= 1
      end
      p[7]==2 ? u[49] += 1 : 0
      p[8]==2 ? u[50] += 1 : 0
  elseif temp<cum_rate[19]
      if maximum([p[7],p[8]])>3
        u[11] -= 1
        u[12] += 1
      else
        u[11] -= 1
      end
      p[7]==3 ? u[49] += 1 : 0
      p[8]==3 ? u[50] += 1 : 0
  elseif temp<cum_rate[20]
      if maximum([p[7],p[8]])>4
        u[12] -= 1
        u[13] += 1
      else
        u[12] -= 1
      end
      p[7]==4 ? u[49] += 1 : 0
      p[8]==4 ? u[50] += 1 : 0
  elseif temp<cum_rate[21]
      if maximum([p[7],p[8]])>5
        u[13] -= 1
        u[14] += 1
      else
        u[13] -= 1
      end
      p[7]==5 ? u[49] += 1 : 0
      p[8]==5 ? u[50] += 1 : 0
  elseif temp<cum_rate[22]
      if maximum([p[7],p[8]])>6
        u[14] -= 1
        u[15] += 1
      else
        u[14] -= 1
      end
      p[7]==6 ? u[49] += 1 : 0
      p[8]==6 ? u[50] += 1 : 0
  elseif temp<cum_rate[23]
      if maximum([p[7],p[8]])>7
        u[15] -= 1
        u[16] += 1
      else
        u[15] -= 1
      end
      p[7]==7 ? u[49] += 1 : 0
      p[8]==7 ? u[50] += 1 : 0
  elseif temp<cum_rate[24]
      if maximum([p[7],p[8]])>8
        u[16] -= 1
        u[17] += 1
      else
        u[16] -= 1
      end
      p[7]==8 ? u[49] += 1 : 0
      p[8]==8 ? u[50] += 1 : 0
  elseif temp<cum_rate[25]
      if maximum([p[7],p[8]])>9
        u[17] -= 1
        u[18] += 1
      else
        u[17] -= 1
      end
      p[7]==9 ? u[49] += 1 : 0
      p[8]==9 ? u[50] += 1 : 0
  elseif temp<cum_rate[26]
      if maximum([p[7],p[8]])>10
        u[18] -= 1
        u[19] += 1
      else
        u[18] -= 1
      end
      p[7]==10 ? u[49] += 1 : 0
      p[8]==10 ? u[50] += 1 : 0
  elseif temp<cum_rate[27]
      if maximum([p[7],p[8]])>11
        u[19] -= 1
        u[20] += 1
      else
        u[19] -= 1
      end
      p[7]==11 ? u[49] += 1 : 0
      p[8]==11 ? u[50] += 1 : 0
  elseif temp<cum_rate[28]
      if maximum([p[7],p[8]])>12
        u[20] -= 1
        u[21] += 1
      else
        u[20] -= 1
      end
      p[7]==12 ? u[49] += 1 : 0
      p[8]==12 ? u[50] += 1 : 0
  elseif temp<cum_rate[29]
      if maximum([p[7],p[8]])>13
        u[21] -= 1
        u[22] += 1
      else
        u[21] -= 1
      end
      p[7]==13 ? u[49] += 1 : 0
      p[8]==13 ? u[50] += 1 : 0
  elseif temp<cum_rate[30]
      if maximum([p[7],p[8]])>14
        u[22] -= 1
        u[23] += 1
      else
        u[22] -= 1
      end
      p[7]==14 ? u[49] += 1 : 0
      p[8]==14 ? u[50] += 1 : 0
  elseif temp<cum_rate[31]
      if maximum([p[7],p[8]])>15
        u[23] -= 1
        u[24] += 1
      else
        u[23] -= 1
      end
      p[7]==15 ? u[49] += 1 : 0
      p[8]==15 ? u[50] += 1 : 0
  elseif temp<cum_rate[32]
      if maximum([p[7],p[8]])>16
        u[24] -= 1
        u[25] += 1
      else
        u[24] -= 1
      end
      p[7]==16 ? u[49] += 1 : 0
      p[8]==16 ? u[50] += 1 : 0
  elseif temp<cum_rate[33]
      if maximum([p[7],p[8]])>17
        u[25] -= 1
        u[26] += 1
      else
        u[25] -= 1
      end
      p[7]==17 ? u[49] += 1 : 0
      p[8]==17 ? u[50] += 1 : 0
  elseif temp<cum_rate[34]
      if maximum([p[7],p[8]])>18
        u[26] -= 1
        u[27] += 1
      else
        u[26] -= 1
      end
      p[7]==18 ? u[49] += 1 : 0
      p[8]==18 ? u[50] += 1 : 0
  elseif temp<cum_rate[35]
      if maximum([p[7],p[8]])>19
        u[27] -= 1
        u[28] += 1
      else
        u[27] -= 1
      end
      p[7]==19 ? u[49] += 1 : 0
      p[8]==19 ? u[50] += 1 : 0
  elseif temp<cum_rate[36]
      u[28] -= 1
      p[7]==20 ? u[49] += 1 : 0
      p[8]==20 ? u[50] += 1 : 0
  elseif temp<cum_rate[37]
      if maximum([p[7],p[8]])>1
        u[30] += 1
        u[2] -= 1
      else
        u[2] -=1
      end
      p[7]==1 ? u[51] += 1 : 0
      p[8]==1 ? u[52] += 1 : 0
  elseif temp<cum_rate[38]
      if maximum([p[7],p[8]])>2
        u[30] -= 1
        u[31] += 1
      else
        u[30] -= 1
      end
      p[7]==2 ? u[51] += 1 : 0
      p[8]==2 ? u[52] += 1 : 0
  elseif temp<cum_rate[39]
      if maximum([p[7],p[8]])>3
        u[31] -= 1
        u[32] += 1
      else
        u[31] -= 1
      end
      p[7]==3 ? u[51] += 1 : 0
      p[8]==3 ? u[52] += 1 : 0
  elseif temp<cum_rate[40]
      if maximum([p[7],p[8]])>4
        u[32] -= 1
        u[33] += 1
      else
        u[32] -= 1
      end
      p[7]==4 ? u[51] += 1 : 0
      p[8]==4 ? u[52] += 1 : 0
  elseif temp<cum_rate[41]
      if maximum([p[7],p[8]])>5
        u[33] -= 1
        u[34] += 1
      else
        u[33] -= 1
      end
      p[7]==5 ? u[51] += 1 : 0
      p[8]==5 ? u[52] += 1 : 0
  elseif temp<cum_rate[42]
      if maximum([p[7],p[8]])>6
        u[34] -= 1
        u[35] += 1
      else
        u[34] -= 1
      end
      p[7]==6 ? u[51] += 1 : 0
      p[8]==6 ? u[52] += 1 : 0
  elseif temp<cum_rate[43]
      if maximum([p[7],p[8]])>7
        u[35] -= 1
        u[36] += 1
      else
        u[35] -= 1
      end
      p[7]==7 ? u[51] += 1 : 0
      p[8]==7 ? u[52] += 1 : 0
  elseif temp<cum_rate[44]
      if maximum([p[7],p[8]])>8
        u[36] -= 1
        u[37] += 1
      else
        u[36] -= 1
      end
      p[7]==8 ? u[51] += 1 : 0
      p[8]==8 ? u[52] += 1 : 0
  elseif temp<cum_rate[45]
      if maximum([p[7],p[8]])>9
        u[37] -= 1
        u[38] += 1
      else
        u[37] -= 1
      end
      p[7]==9 ? u[51] += 1 : 0
      p[8]==9 ? u[52] += 1 : 0
  elseif temp<cum_rate[46]
      if maximum([p[7],p[8]])>10
        u[38] -= 1
        u[39] += 1
      else
        u[38] -= 1
      end
      p[7]==10 ? u[51] += 1 : 0
      p[8]==10 ? u[52] += 1 : 0
  elseif temp<cum_rate[47]
      if maximum([p[7],p[8]])>11
        u[39] -= 1
        u[40] += 1
      else
        u[39] -= 1
      end
      p[7]==11 ? u[51] += 1 : 0
      p[8]==11 ? u[52] += 1 : 0
  elseif temp<cum_rate[48]
      if maximum([p[7],p[8]])>12
        u[40] -= 1
        u[41] += 1
      else
        u[40] -= 1
      end
      p[7]==12 ? u[51] += 1 : 0
      p[8]==12 ? u[52] += 1 : 0
  elseif temp<cum_rate[49]
      if maximum([p[7],p[8]])>13
        u[41] -= 1
        u[42] += 1
      else
        u[41] -= 1
      end
      p[7]==13 ? u[51] += 1 : 0
      p[8]==13 ? u[52] += 1 : 0
  elseif temp<cum_rate[50]
      if maximum([p[7],p[8]])>14
        u[42] -= 1
        u[43] += 1
      else
        u[42] -= 1
      end
      p[7]==14 ? u[51] += 1 : 0
      p[8]==14 ? u[52] += 1 : 0
  elseif temp<cum_rate[51]
      if maximum([p[7],p[8]])>15
        u[43] -= 1
        u[44] += 1
      else
        u[43] -= 1
      end
      p[7]==15 ? u[51] += 1 : 0
      p[8]==15 ? u[52] += 1 : 0
  elseif temp<cum_rate[52]
      if maximum([p[7],p[8]])>16
        u[44] -= 1
        u[45] += 1
      else
        u[44] -= 1
      end
      p[7]==16 ? u[51] += 1 : 0
      p[8]==16 ? u[52] += 1 : 0
  elseif temp<cum_rate[53]
      if maximum([p[7],p[8]])>17
        u[45] -= 1
        u[46] += 1
      else
        u[45] -= 1
      end
      p[7]==17 ? u[51] += 1 : 0
      p[8]==17 ? u[52] += 1 : 0
  elseif temp<cum_rate[54]
      if maximum([p[7],p[8]])>18
        u[46] -= 1
        u[47] += 1
      else
        u[46] -= 1
      end
      p[7]==18 ? u[51] += 1 : 0
      p[8]==18 ? u[52] += 1 : 0
  elseif temp<cum_rate[55]
      if maximum([p[7],p[8]])>19
        u[47] -= 1
        u[48] += 1
      else
        u[47] -= 1
      end
      p[7]==19 ? u[51] += 1 : 0
      p[8]==19 ? u[52] += 1 : 0
  elseif temp<cum_rate[56]
      u[48] -= 1
      p[7]==20 ? u[51] += 1 : 0
      p[8]==20 ? u[52] += 1 : 0
  end
    t = t+t_step
    push!(out,copy(u))
    push!(t_out,copy(t))
  end
  return out, t_out
end

### Dox induction experiment
function sim_delay_dox(p::Array{Float64,1},tspan::Tuple{Float64,Float64},u0::Array{Int64,1})
  t::Float64=tspan[1]
  i::Int64=0
  out::Array{Array{Int64,1}} = []
  t_out::Array{Float64,1} = []
  push!(out, u0)
  push!(t_out,t)
  u::Array{Int64,1}=copy(u0)
  #############################################################
  rate::Array{Float64,1}=zeros(56)
  while (t<tspan[2])
    i+=1
    if t<24
        p[24] = 10
        p[29] = 0
        p[25] = 0
        p[30] = 0
    else
        p[24] = 10
        p[29] = 0
        p[25] = 0
        p[30] = 1
    end

    # Xist transcription
    rate[1] = p[21]*(p[24]+p[29]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[7]^p[13]/(u[7]^p[13]+(p[22]*p[14])^p[13])))
    rate[2] = p[21]*(p[25]+p[30]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[8]^p[13]/(u[8]^p[13]+(p[22]*p[14])^p[13])))
    # cXR production
    rate[3] = p[22]*(1-u[51]^p[5]/(u[51]^p[5]+(p[21]*p[6])^p[5]))
    rate[4] = p[22]*(1-u[49]^p[5]/(u[49]^p[5]+(p[21]*p[6])^p[5]))
    # tXA production
    rate[5] = p[23]*(1-u[52]^p[3]/(u[52]^p[3]+(p[21]*p[4])^p[3]))
    rate[6] = p[23]*(1-u[50]^p[3]/(u[50]^p[3]+(p[21]*p[4])^p[3]))
    # degradation
    rate[7] = 0.1733*u[3]
    rate[8] = 0.1733*u[4]
    rate[9] = 1*u[7]
    rate[10] = 1*u[8]
    rate[11] = 1*u[5]
    rate[12] = 1*u[6]
    rate[13] = 1*u[49]
    rate[14] = 1*u[50]
    rate[15] = 1*u[51]
    rate[16] = 1*u[52]
    ### silencing intermediates on chr1
    rate[17] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[1],0)
    rate[18] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[10],0)
    rate[19] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[11],0)
    rate[20] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[12],0)
    rate[21] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[13],0)
    rate[22] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[14],0)
    rate[23] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[15],0)
    rate[24] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[16],0)
    rate[25] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[17],0)
    rate[26] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[18],0)
    rate[27] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[19],0)
    rate[28] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[20],0)
    rate[29] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[21],0)
    rate[30] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[22],0)
    rate[31] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[23],0)
    rate[32] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[24],0)
    rate[33] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[25],0)
    rate[34] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[26],0)
    rate[35] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[27],0)
    rate[36] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[28],0)
    ### silencing intermediates chr 2
    rate[37] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[2],0)
    rate[38] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[30],0)
    rate[39] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[31],0)
    rate[40] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[32],0)
    rate[41] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[33],0)
    rate[42] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[34],0)
    rate[43] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[35],0)
    rate[44] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[36],0)
    rate[45] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[37],0)
    rate[46] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[38],0)
    rate[47] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[39],0)
    rate[48] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[40],0)
    rate[49] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[41],0)
    rate[50] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[42],0)
    rate[51] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[43],0)
    rate[52] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[44],0)
    rate[53] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[45],0)
    rate[54] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[46],0)
    rate[55] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[47],0)
    rate[56] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[48],0)

    # Gillespie
    r = rand(2)
    # determine time point of next reaction
    t_step::Float64 = 1/sum(rate)*log(1/(1-r[1]))
    # decide which reaction will occur
    temp::Float64 = r[2]*sum(rate)
    cum_rate = cumsum(rate)
    # Execute reaction
    if temp<cum_rate[1]
      u[1] += 1
      u[3] += 1
        if p[7]==0
        u[49] += 1
      end
      if p[8]==0
        u[50] += 1
      end
  elseif temp<cum_rate[2]
      u[2] += 1
      u[4] += 1
      if p[7]==0
        u[51] += 1
      end
      if p[8]==0
        u[52] += 1
      end
  elseif temp<cum_rate[3]
        u[8] += 1
  elseif temp<cum_rate[4]
      u[7] += 1
  elseif temp<cum_rate[5]
        u[6] += 1
  elseif temp<cum_rate[6]
        u[5] += 1
  elseif temp<cum_rate[7]
      if maximum([p[7],p[8]])==0
        u[1] -= 1
      end
        u[3] -=1
  elseif temp<cum_rate[8]
      if maximum([p[7],p[8]])==0
        u[2] -= 1
      end
        u[4] -=1
  elseif temp<cum_rate[9]
        u[7] -= 1
  elseif temp<cum_rate[10]
        u[8] -= 1
  elseif temp<cum_rate[11]
        u[5] -= 1
  elseif temp<cum_rate[12]
        u[6] -= 1
  elseif temp<cum_rate[13]
        u[49] -= 1
  elseif temp<cum_rate[14]
        u[50] -= 1
  elseif temp<cum_rate[15]
        u[51] -= 1
  elseif temp<cum_rate[16]
        u[52] -= 1
  elseif temp<cum_rate[17]
      if maximum([p[7],p[8]])>1
        u[10] += 1
        u[1] -= 1
      else
        u[1] -= 1
      end
      p[7]==1 ? u[49] += 1 : 0
      p[8]==1 ? u[50] += 1 : 0
  elseif temp<cum_rate[18]
      if maximum([p[7],p[8]])>2
        u[10] -= 1
        u[11] += 1
      else
        u[10] -= 1
      end
      p[7]==2 ? u[49] += 1 : 0
      p[8]==2 ? u[50] += 1 : 0
  elseif temp<cum_rate[19]
      if maximum([p[7],p[8]])>3
        u[11] -= 1
        u[12] += 1
      else
        u[11] -= 1
      end
      p[7]==3 ? u[49] += 1 : 0
      p[8]==3 ? u[50] += 1 : 0
  elseif temp<cum_rate[20]
      if maximum([p[7],p[8]])>4
        u[12] -= 1
        u[13] += 1
      else
        u[12] -= 1
      end
      p[7]==4 ? u[49] += 1 : 0
      p[8]==4 ? u[50] += 1 : 0
  elseif temp<cum_rate[21]
      if maximum([p[7],p[8]])>5
        u[13] -= 1
        u[14] += 1
      else
        u[13] -= 1
      end
      p[7]==5 ? u[49] += 1 : 0
      p[8]==5 ? u[50] += 1 : 0
  elseif temp<cum_rate[22]
      if maximum([p[7],p[8]])>6
        u[14] -= 1
        u[15] += 1
      else
        u[14] -= 1
      end
      p[7]==6 ? u[49] += 1 : 0
      p[8]==6 ? u[50] += 1 : 0
  elseif temp<cum_rate[23]
      if maximum([p[7],p[8]])>7
        u[15] -= 1
        u[16] += 1
      else
        u[15] -= 1
      end
      p[7]==7 ? u[49] += 1 : 0
      p[8]==7 ? u[50] += 1 : 0
  elseif temp<cum_rate[24]
      if maximum([p[7],p[8]])>8
        u[16] -= 1
        u[17] += 1
      else
        u[16] -= 1
      end
      p[7]==8 ? u[49] += 1 : 0
      p[8]==8 ? u[50] += 1 : 0
  elseif temp<cum_rate[25]
      if maximum([p[7],p[8]])>9
        u[17] -= 1
        u[18] += 1
      else
        u[17] -= 1
      end
      p[7]==9 ? u[49] += 1 : 0
      p[8]==9 ? u[50] += 1 : 0
  elseif temp<cum_rate[26]
      if maximum([p[7],p[8]])>10
        u[18] -= 1
        u[19] += 1
      else
        u[18] -= 1
      end
      p[7]==10 ? u[49] += 1 : 0
      p[8]==10 ? u[50] += 1 : 0
  elseif temp<cum_rate[27]
      if maximum([p[7],p[8]])>11
        u[19] -= 1
        u[20] += 1
      else
        u[19] -= 1
      end
      p[7]==11 ? u[49] += 1 : 0
      p[8]==11 ? u[50] += 1 : 0
  elseif temp<cum_rate[28]
      if maximum([p[7],p[8]])>12
        u[20] -= 1
        u[21] += 1
      else
        u[20] -= 1
      end
      p[7]==12 ? u[49] += 1 : 0
      p[8]==12 ? u[50] += 1 : 0
  elseif temp<cum_rate[29]
      if maximum([p[7],p[8]])>13
        u[21] -= 1
        u[22] += 1
      else
        u[21] -= 1
      end
      p[7]==13 ? u[49] += 1 : 0
      p[8]==13 ? u[50] += 1 : 0
  elseif temp<cum_rate[30]
      if maximum([p[7],p[8]])>14
        u[22] -= 1
        u[23] += 1
      else
        u[22] -= 1
      end
      p[7]==14 ? u[49] += 1 : 0
      p[8]==14 ? u[50] += 1 : 0
  elseif temp<cum_rate[31]
      if maximum([p[7],p[8]])>15
        u[23] -= 1
        u[24] += 1
      else
        u[23] -= 1
      end
      p[7]==15 ? u[49] += 1 : 0
      p[8]==15 ? u[50] += 1 : 0
  elseif temp<cum_rate[32]
      if maximum([p[7],p[8]])>16
        u[24] -= 1
        u[25] += 1
      else
        u[24] -= 1
      end
      p[7]==16 ? u[49] += 1 : 0
      p[8]==16 ? u[50] += 1 : 0
  elseif temp<cum_rate[33]
      if maximum([p[7],p[8]])>17
        u[25] -= 1
        u[26] += 1
      else
        u[25] -= 1
      end
      p[7]==17 ? u[49] += 1 : 0
      p[8]==17 ? u[50] += 1 : 0
  elseif temp<cum_rate[34]
      if maximum([p[7],p[8]])>18
        u[26] -= 1
        u[27] += 1
      else
        u[26] -= 1
      end
      p[7]==18 ? u[49] += 1 : 0
      p[8]==18 ? u[50] += 1 : 0
  elseif temp<cum_rate[35]
      if maximum([p[7],p[8]])>19
        u[27] -= 1
        u[28] += 1
      else
        u[27] -= 1
      end
      p[7]==19 ? u[49] += 1 : 0
      p[8]==19 ? u[50] += 1 : 0
  elseif temp<cum_rate[36]
      u[28] -= 1
      p[7]==20 ? u[49] += 1 : 0
      p[8]==20 ? u[50] += 1 : 0
  elseif temp<cum_rate[37]
      if maximum([p[7],p[8]])>1
        u[30] += 1
        u[2] -= 1
      else
        u[2] -= 1
      end
      p[7]==1 ? u[51] += 1 : 0
      p[8]==1 ? u[52] += 1 : 0
  elseif temp<cum_rate[38]
      if maximum([p[7],p[8]])>2
        u[30] -= 1
        u[31] += 1
      else
        u[30] -= 1
      end
      p[7]==2 ? u[51] += 1 : 0
      p[8]==2 ? u[52] += 1 : 0
  elseif temp<cum_rate[39]
      if maximum([p[7],p[8]])>3
        u[31] -= 1
        u[32] += 1
      else
        u[31] -= 1
      end
      p[7]==3 ? u[51] += 1 : 0
      p[8]==3 ? u[52] += 1 : 0
  elseif temp<cum_rate[40]
      if maximum([p[7],p[8]])>4
        u[32] -= 1
        u[33] += 1
      else
        u[32] -= 1
      end
      p[7]==4 ? u[51] += 1 : 0
      p[8]==4 ? u[52] += 1 : 0
  elseif temp<cum_rate[41]
      if maximum([p[7],p[8]])>5
        u[33] -= 1
        u[34] += 1
      else
        u[33] -= 1
      end
      p[7]==5 ? u[51] += 1 : 0
      p[8]==5 ? u[52] += 1 : 0
  elseif temp<cum_rate[42]
      if maximum([p[7],p[8]])>6
        u[34] -= 1
        u[35] += 1
      else
        u[34] -= 1
      end
      p[7]==6 ? u[51] += 1 : 0
      p[8]==6 ? u[52] += 1 : 0
  elseif temp<cum_rate[43]
      if maximum([p[7],p[8]])>7
        u[35] -= 1
        u[36] += 1
      else
        u[35] -= 1
      end
      p[7]==7 ? u[51] += 1 : 0
      p[8]==7 ? u[52] += 1 : 0
  elseif temp<cum_rate[44]
      if maximum([p[7],p[8]])>8
        u[36] -= 1
        u[37] += 1
      else
        u[36] -= 1
      end
      p[7]==8 ? u[51] += 1 : 0
      p[8]==8 ? u[52] += 1 : 0
  elseif temp<cum_rate[45]
      if maximum([p[7],p[8]])>9
        u[37] -= 1
        u[38] += 1
      else
        u[37] -= 1
      end
      p[7]==9 ? u[51] += 1 : 0
      p[8]==9 ? u[52] += 1 : 0
  elseif temp<cum_rate[46]
      if maximum([p[7],p[8]])>10
        u[38] -= 1
        u[39] += 1
      else
        u[38] -= 1
      end
      p[7]==10 ? u[51] += 1 : 0
      p[8]==10 ? u[52] += 1 : 0
  elseif temp<cum_rate[47]
      if maximum([p[7],p[8]])>11
        u[39] -= 1
        u[40] += 1
      else
        u[39] -= 1
      end
      p[7]==11 ? u[51] += 1 : 0
      p[8]==11 ? u[52] += 1 : 0
  elseif temp<cum_rate[48]
      if maximum([p[7],p[8]])>12
        u[40] -= 1
        u[41] += 1
      else
        u[40] -= 1
      end
      p[7]==12 ? u[51] += 1 : 0
      p[8]==12 ? u[52] += 1 : 0
  elseif temp<cum_rate[49]
      if maximum([p[7],p[8]])>13
        u[41] -= 1
        u[42] += 1
      else
        u[41] -= 1
      end
      p[7]==13 ? u[51] += 1 : 0
      p[8]==13 ? u[52] += 1 : 0
  elseif temp<cum_rate[50]
      if maximum([p[7],p[8]])>14
        u[42] -= 1
        u[43] += 1
      else
        u[42] -= 1
      end
      p[7]==14 ? u[51] += 1 : 0
      p[8]==14 ? u[52] += 1 : 0
  elseif temp<cum_rate[51]
      if maximum([p[7],p[8]])>15
        u[43] -= 1
        u[44] += 1
      else
        u[43] -= 1
      end
      p[7]==15 ? u[51] += 1 : 0
      p[8]==15 ? u[52] += 1 : 0
  elseif temp<cum_rate[52]
      if maximum([p[7],p[8]])>16
        u[44] -= 1
        u[45] += 1
      else
        u[44] -= 1
      end
      p[7]==16 ? u[51] += 1 : 0
      p[8]==16 ? u[52] += 1 : 0
  elseif temp<cum_rate[53]
      if maximum([p[7],p[8]])>17
        u[45] -= 1
        u[46] += 1
      else
        u[45] -= 1
      end
      p[7]==17 ? u[51] += 1 : 0
      p[8]==17 ? u[52] += 1 : 0
  elseif temp<cum_rate[54]
      if maximum([p[7],p[8]])>18
        u[46] -= 1
        u[47] += 1
      else
        u[46] -= 1
      end
      p[7]==18 ? u[51] += 1 : 0
      p[8]==18 ? u[52] += 1 : 0
  elseif temp<cum_rate[55]
      if maximum([p[7],p[8]])>19
        u[47] -= 1
        u[48] += 1
      else
        u[47] -= 1
      end
      p[7]==19 ? u[51] += 1 : 0
      p[8]==19 ? u[52] += 1 : 0
  elseif temp<cum_rate[56]
      u[48] -= 1
      p[7]==20 ? u[51] += 1 : 0
      p[8]==20 ? u[52] += 1 : 0
  end
    t = t+t_step
    push!(out,copy(u))
    push!(t_out,copy(t))
  end
  return out, t_out
end

### Artificial biallelic induction experiment
function sim_delay_ba(p::Array{Float64,1},tspan::Tuple{Float64,Float64},u0::Array{Int64,1})
  t::Float64=tspan[1]
  i::Int64=0
  out::Array{Array{Int64,1}} = []
  t_out::Array{Float64,1} = []
  push!(out, u0)
  push!(t_out,t)
  u::Array{Int64,1}=copy(u0)
  #############################################################
  rate::Array{Float64,1}=zeros(56)
  while (t<tspan[2])
    i+=1
    if t<48
        p[24] = 0
        p[29] = 1/3
        p[25] = 0
        p[30] = 1
    else
        p[24] = 10
        p[29] = 0
        p[25] = 0
        p[30] = 1
    end

    # Xist transcription
    rate[1] = p[21]*(p[24]+p[29]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[7]^p[13]/(u[7]^p[13]+(p[22]*p[14])^p[13])))
    rate[2] = p[21]*(p[25]+p[30]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[8]^p[13]/(u[8]^p[13]+(p[22]*p[14])^p[13])))
    # cXR production
    rate[3] = p[22]*(1-u[51]^p[5]/(u[51]^p[5]+(p[21]*p[6])^p[5]))
    rate[4] = p[22]*(1-u[49]^p[5]/(u[49]^p[5]+(p[21]*p[6])^p[5]))
    # tXA production
    rate[5] = p[23]*(1-u[52]^p[3]/(u[52]^p[3]+(p[21]*p[4])^p[3]))
    rate[6] = p[23]*(1-u[50]^p[3]/(u[50]^p[3]+(p[21]*p[4])^p[3]))
    # degradation
    rate[7] = 0.1733*u[3]
    rate[8] = 0.1733*u[4]
    rate[9] = 1*u[7]
    rate[10] = 1*u[8]
    rate[11] = 1*u[5]
    rate[12] = 1*u[6]
    rate[13] = 1*u[49]
    rate[14] = 1*u[50]
    rate[15] = 1*u[51]
    rate[16] = 1*u[52]
    ### silencing intermediates on chr1
    rate[17] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[1],0)
    rate[18] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[10],0)
    rate[19] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[11],0)
    rate[20] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[12],0)
    rate[21] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[13],0)
    rate[22] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[14],0)
    rate[23] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[15],0)
    rate[24] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[16],0)
    rate[25] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[17],0)
    rate[26] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[18],0)
    rate[27] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[19],0)
    rate[28] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[20],0)
    rate[29] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[21],0)
    rate[30] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[22],0)
    rate[31] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[23],0)
    rate[32] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[24],0)
    rate[33] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[25],0)
    rate[34] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[26],0)
    rate[35] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[27],0)
    rate[36] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[28],0)
    ### silencing intermediates chr 2
    rate[37] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[2],0)
    rate[38] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[30],0)
    rate[39] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[31],0)
    rate[40] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[32],0)
    rate[41] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[33],0)
    rate[42] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[34],0)
    rate[43] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[35],0)
    rate[44] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[36],0)
    rate[45] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[37],0)
    rate[46] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[38],0)
    rate[47] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[39],0)
    rate[48] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[40],0)
    rate[49] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[41],0)
    rate[50] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[42],0)
    rate[51] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[43],0)
    rate[52] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[44],0)
    rate[53] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[45],0)
    rate[54] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[46],0)
    rate[55] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[47],0)
    rate[56] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[48],0)

	# Gillespie
    r = rand(2)
    # determine time point of next reaction
    t_step::Float64 = 1/sum(rate)*log(1/(1-r[1]))
    # decide which reaction will occur
    temp::Float64 = r[2]*sum(rate)
    cum_rate = cumsum(rate)
    # Execute reaction
    if temp<cum_rate[1]
      u[1] += 1
      u[3] += 1
        if p[7]==0
           u[49] += 1
        end
        if p[8]==0
           u[50] += 1
        end
  elseif temp<cum_rate[2]
      u[2] += 1
      u[4] += 1
      if p[7]==0
        u[51] += 1
      end
      if p[8]==0
        u[52] += 1
      end
  elseif temp<cum_rate[3]
        u[8] += 1
  elseif temp<cum_rate[4]
      u[7] += 1
  elseif temp<cum_rate[5]
        u[6] += 1
  elseif temp<cum_rate[6]
        u[5] += 1
  elseif temp<cum_rate[7]
      if maximum([p[7],p[8]])==0
        u[1] -= 1
      end
        u[3] -=1
  elseif temp<cum_rate[8]
      if maximum([p[7],p[8]])==0
        u[2] -= 1
      end
        u[4] -=1
  elseif temp<cum_rate[9]
        u[7] -= 1
  elseif temp<cum_rate[10]
        u[8] -= 1
  elseif temp<cum_rate[11]
        u[5] -= 1
  elseif temp<cum_rate[12]
        u[6] -= 1
  elseif temp<cum_rate[13]
        u[49] -= 1
  elseif temp<cum_rate[14]
        u[50] -= 1
  elseif temp<cum_rate[15]
        u[51] -= 1
  elseif temp<cum_rate[16]
        u[52] -= 1
  elseif temp<cum_rate[17]
      if maximum([p[7],p[8]])>1
        u[10] += 1
        u[1] -= 1
      else
        u[1] -= 1
      end
      p[7]==1 ? u[49] += 1 : 0
      p[8]==1 ? u[50] += 1 : 0
  elseif temp<cum_rate[18]
      if maximum([p[7],p[8]])>2
        u[10] -= 1
        u[11] += 1
      else
        u[10] -= 1
      end
      p[7]==2 ? u[49] += 1 : 0
      p[8]==2 ? u[50] += 1 : 0
  elseif temp<cum_rate[19]
      if maximum([p[7],p[8]])>3
        u[11] -= 1
        u[12] += 1
      else
        u[11] -= 1
      end
      p[7]==3 ? u[49] += 1 : 0
      p[8]==3 ? u[50] += 1 : 0
  elseif temp<cum_rate[20]
      if maximum([p[7],p[8]])>4
        u[12] -= 1
        u[13] += 1
      else
        u[12] -= 1
      end
      p[7]==4 ? u[49] += 1 : 0
      p[8]==4 ? u[50] += 1 : 0
  elseif temp<cum_rate[21]
      if maximum([p[7],p[8]])>5
        u[13] -= 1
        u[14] += 1
      else
        u[13] -= 1
      end
      p[7]==5 ? u[49] += 1 : 0
      p[8]==5 ? u[50] += 1 : 0
  elseif temp<cum_rate[22]
      if maximum([p[7],p[8]])>6
        u[14] -= 1
        u[15] += 1
      else
        u[14] -= 1
      end
      p[7]==6 ? u[49] += 1 : 0
      p[8]==6 ? u[50] += 1 : 0
  elseif temp<cum_rate[23]
      if maximum([p[7],p[8]])>7
        u[15] -= 1
        u[16] += 1
      else
        u[15] -= 1
      end
      p[7]==7 ? u[49] += 1 : 0
      p[8]==7 ? u[50] += 1 : 0
  elseif temp<cum_rate[24]
      if maximum([p[7],p[8]])>8
        u[16] -= 1
        u[17] += 1
      else
        u[16] -= 1
      end
      p[7]==8 ? u[49] += 1 : 0
      p[8]==8 ? u[50] += 1 : 0
  elseif temp<cum_rate[25]
      if maximum([p[7],p[8]])>9
        u[17] -= 1
        u[18] += 1
      else
        u[17] -= 1
      end
      p[7]==9 ? u[49] += 1 : 0
      p[8]==9 ? u[50] += 1 : 0
  elseif temp<cum_rate[26]
      if maximum([p[7],p[8]])>10
        u[18] -= 1
        u[19] += 1
      else
        u[18] -= 1
      end
      p[7]==10 ? u[49] += 1 : 0
      p[8]==10 ? u[50] += 1 : 0
  elseif temp<cum_rate[27]
      if maximum([p[7],p[8]])>11
        u[19] -= 1
        u[20] += 1
      else
        u[19] -= 1
      end
      p[7]==11 ? u[49] += 1 : 0
      p[8]==11 ? u[50] += 1 : 0
  elseif temp<cum_rate[28]
      if maximum([p[7],p[8]])>12
        u[20] -= 1
        u[21] += 1
      else
        u[20] -= 1
      end
      p[7]==12 ? u[49] += 1 : 0
      p[8]==12 ? u[50] += 1 : 0
  elseif temp<cum_rate[29]
      if maximum([p[7],p[8]])>13
        u[21] -= 1
        u[22] += 1
      else
        u[21] -= 1
      end
      p[7]==13 ? u[49] += 1 : 0
      p[8]==13 ? u[50] += 1 : 0
  elseif temp<cum_rate[30]
      if maximum([p[7],p[8]])>14
        u[22] -= 1
        u[23] += 1
      else
        u[22] -= 1
      end
      p[7]==14 ? u[49] += 1 : 0
      p[8]==14 ? u[50] += 1 : 0
  elseif temp<cum_rate[31]
      if maximum([p[7],p[8]])>15
        u[23] -= 1
        u[24] += 1
      else
        u[23] -= 1
      end
      p[7]==15 ? u[49] += 1 : 0
      p[8]==15 ? u[50] += 1 : 0
  elseif temp<cum_rate[32]
      if maximum([p[7],p[8]])>16
        u[24] -= 1
        u[25] += 1
      else
        u[24] -= 1
      end
      p[7]==16 ? u[49] += 1 : 0
      p[8]==16 ? u[50] += 1 : 0
  elseif temp<cum_rate[33]
      if maximum([p[7],p[8]])>17
        u[25] -= 1
        u[26] += 1
      else
        u[25] -= 1
      end
      p[7]==17 ? u[49] += 1 : 0
      p[8]==17 ? u[50] += 1 : 0
  elseif temp<cum_rate[34]
      if maximum([p[7],p[8]])>18
        u[26] -= 1
        u[27] += 1
      else
        u[26] -= 1
      end
      p[7]==18 ? u[49] += 1 : 0
      p[8]==18 ? u[50] += 1 : 0
  elseif temp<cum_rate[35]
      if maximum([p[7],p[8]])>19
        u[27] -= 1
        u[28] += 1
      else
        u[27] -= 1
      end
      p[7]==19 ? u[49] += 1 : 0
      p[8]==19 ? u[50] += 1 : 0
  elseif temp<cum_rate[36]
      u[28] -= 1
      p[7]==20 ? u[49] += 1 : 0
      p[8]==20 ? u[50] += 1 : 0
  elseif temp<cum_rate[37]
      if maximum([p[7],p[8]])>1
        u[30] += 1
        u[2] -= 1
      else
        u[2] -=1
      end
      p[7]==1 ? u[51] += 1 : 0
      p[8]==1 ? u[52] += 1 : 0
  elseif temp<cum_rate[38]
      if maximum([p[7],p[8]])>2
        u[30] -= 1
        u[31] += 1
      else
        u[30] -= 1
      end
      p[7]==2 ? u[51] += 1 : 0
      p[8]==2 ? u[52] += 1 : 0
  elseif temp<cum_rate[39]
      if maximum([p[7],p[8]])>3
        u[31] -= 1
        u[32] += 1
      else
        u[31] -= 1
      end
      p[7]==3 ? u[51] += 1 : 0
      p[8]==3 ? u[52] += 1 : 0
  elseif temp<cum_rate[40]
      if maximum([p[7],p[8]])>4
        u[32] -= 1
        u[33] += 1
      else
        u[32] -= 1
      end
      p[7]==4 ? u[51] += 1 : 0
      p[8]==4 ? u[52] += 1 : 0
  elseif temp<cum_rate[41]
      if maximum([p[7],p[8]])>5
        u[33] -= 1
        u[34] += 1
      else
        u[33] -= 1
      end
      p[7]==5 ? u[51] += 1 : 0
      p[8]==5 ? u[52] += 1 : 0
  elseif temp<cum_rate[42]
      if maximum([p[7],p[8]])>6
        u[34] -= 1
        u[35] += 1
      else
        u[34] -= 1
      end
      p[7]==6 ? u[51] += 1 : 0
      p[8]==6 ? u[52] += 1 : 0
  elseif temp<cum_rate[43]
      if maximum([p[7],p[8]])>7
        u[35] -= 1
        u[36] += 1
      else
        u[35] -= 1
      end
      p[7]==7 ? u[51] += 1 : 0
      p[8]==7 ? u[52] += 1 : 0
  elseif temp<cum_rate[44]
      if maximum([p[7],p[8]])>8
        u[36] -= 1
        u[37] += 1
      else
        u[36] -= 1
      end
      p[7]==8 ? u[51] += 1 : 0
      p[8]==8 ? u[52] += 1 : 0
  elseif temp<cum_rate[45]
      if maximum([p[7],p[8]])>9
        u[37] -= 1
        u[38] += 1
      else
        u[37] -= 1
      end
      p[7]==9 ? u[51] += 1 : 0
      p[8]==9 ? u[52] += 1 : 0
  elseif temp<cum_rate[46]
      if maximum([p[7],p[8]])>10
        u[38] -= 1
        u[39] += 1
      else
        u[38] -= 1
      end
      p[7]==10 ? u[51] += 1 : 0
      p[8]==10 ? u[52] += 1 : 0
  elseif temp<cum_rate[47]
      if maximum([p[7],p[8]])>11
        u[39] -= 1
        u[40] += 1
      else
        u[39] -= 1
      end
      p[7]==11 ? u[51] += 1 : 0
      p[8]==11 ? u[52] += 1 : 0
  elseif temp<cum_rate[48]
      if maximum([p[7],p[8]])>12
        u[40] -= 1
        u[41] += 1
      else
        u[40] -= 1
      end
      p[7]==12 ? u[51] += 1 : 0
      p[8]==12 ? u[52] += 1 : 0
  elseif temp<cum_rate[49]
      if maximum([p[7],p[8]])>13
        u[41] -= 1
        u[42] += 1
      else
        u[41] -= 1
      end
      p[7]==13 ? u[51] += 1 : 0
      p[8]==13 ? u[52] += 1 : 0
  elseif temp<cum_rate[50]
      if maximum([p[7],p[8]])>14
        u[42] -= 1
        u[43] += 1
      else
        u[42] -= 1
      end
      p[7]==14 ? u[51] += 1 : 0
      p[8]==14 ? u[52] += 1 : 0
  elseif temp<cum_rate[51]
      if maximum([p[7],p[8]])>15
        u[43] -= 1
        u[44] += 1
      else
        u[43] -= 1
      end
      p[7]==15 ? u[51] += 1 : 0
      p[8]==15 ? u[52] += 1 : 0
  elseif temp<cum_rate[52]
      if maximum([p[7],p[8]])>16
        u[44] -= 1
        u[45] += 1
      else
        u[44] -= 1
      end
      p[7]==16 ? u[51] += 1 : 0
      p[8]==16 ? u[52] += 1 : 0
  elseif temp<cum_rate[53]
      if maximum([p[7],p[8]])>17
        u[45] -= 1
        u[46] += 1
      else
        u[45] -= 1
      end
      p[7]==17 ? u[51] += 1 : 0
      p[8]==17 ? u[52] += 1 : 0
  elseif temp<cum_rate[54]
      if maximum([p[7],p[8]])>18
        u[46] -= 1
        u[47] += 1
      else
        u[46] -= 1
      end
      p[7]==18 ? u[51] += 1 : 0
      p[8]==18 ? u[52] += 1 : 0
  elseif temp<cum_rate[55]
      if maximum([p[7],p[8]])>19
        u[47] -= 1
        u[48] += 1
      else
        u[47] -= 1
      end
      p[7]==19 ? u[51] += 1 : 0
      p[8]==19 ? u[52] += 1 : 0
  elseif temp<cum_rate[56]
      u[48] -= 1
      p[7]==20 ? u[51] += 1 : 0
      p[8]==20 ? u[52] += 1 : 0
  end
    t = t+t_step
    push!(out,copy(u))
    push!(t_out,copy(t))
  end
  return out, t_out
end

### Heterozygous Tsix deletion
function sim_delay_skew(p::Array{Float64,1},tspan::Tuple{Float64,Float64},u0::Array{Int64,1})
  t::Float64=tspan[1]
  i::Int64=0
  out::Array{Array{Int64,1}} = []
  t_out::Array{Float64,1} = []
  push!(out, u0)
  push!(t_out,t)
  u::Array{Int64,1}=copy(u0)
  #############################################################
  rate::Array{Float64,1}=zeros(56)
  while (t<tspan[2])
    i+=1
    # Xist transcription
    rate[1] = p[21]*(p[24]+p[29]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[7]^p[13]/(u[7]^p[13]+(p[22]*p[14])^p[13])))
    rate[2] = p[21]*(p[25]+p[30]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[8]^p[13]/(u[8]^p[13]+(p[22]*p[14])^p[13])))
    # cXR production
    rate[3] = p[22]*(1-u[51]^p[5]/(u[51]^p[5]+(p[21]*p[6])^p[5]))
    rate[4] = p[22]*(1-u[49]^p[5]/(u[49]^p[5]+(p[21]*p[6])^p[5]))
    # tXA production
    rate[5] = p[23]*(1-u[52]^p[3]/(u[52]^p[3]+(p[21]*p[4])^p[3]))
    rate[6] = p[23]*(1-u[50]^p[3]/(u[50]^p[3]+(p[21]*p[4])^p[3]))
    # degradation
    rate[7] = 0.1733*u[3]
    rate[8] = 0.1733*u[4]
    rate[9] = 1*u[7]
    rate[10] = 1*u[8]
    rate[11] = 1*u[5]
    rate[12] = 1*u[6]
    rate[13] = 1*u[49]
    rate[14] = 1*u[50]
    rate[15] = 1*u[51]
    rate[16] = 1*u[52]
    ### silencing intermediates on chr1
    rate[17] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[1],0)
    rate[18] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[10],0)
    rate[19] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[11],0)
    rate[20] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[12],0)
    rate[21] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[13],0)
    rate[22] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[14],0)
    rate[23] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[15],0)
    rate[24] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[16],0)
    rate[25] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[17],0)
    rate[26] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[18],0)
    rate[27] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[19],0)
    rate[28] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[20],0)
    rate[29] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[21],0)
    rate[30] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[22],0)
    rate[31] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[23],0)
    rate[32] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[24],0)
    rate[33] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[25],0)
    rate[34] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[26],0)
    rate[35] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[27],0)
    rate[36] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[28],0)
    ### silencing intermediates chr 2
    rate[37] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[2],0)
    rate[38] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[30],0)
    rate[39] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[31],0)
    rate[40] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[32],0)
    rate[41] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[33],0)
    rate[42] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[34],0)
    rate[43] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[35],0)
    rate[44] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[36],0)
    rate[45] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[37],0)
    rate[46] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[38],0)
    rate[47] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[39],0)
    rate[48] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[40],0)
    rate[49] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[41],0)
    rate[50] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[42],0)
    rate[51] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[43],0)
    rate[52] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[44],0)
    rate[53] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[45],0)
    rate[54] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[46],0)
    rate[55] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[47],0)
    rate[56] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[48],0)

	# Gillespie
    r = rand(2)
    # determine time point of next reaction
    t_step::Float64 = 1/sum(rate)*log(1/(1-r[1]))
    # decide which reaction will occur
    temp::Float64 = r[2]*sum(rate)
    cum_rate = cumsum(rate)
    # Execute reaction
    if temp<cum_rate[1]
      u[1] += 1
      u[3] += 1
        if p[7]==0
           u[49] += 1
        end
        if p[8]==0
          u[50] += 1
        end
  elseif temp<cum_rate[2]
      u[2] += 1
      u[4] += 1
      if p[7]==0
        u[51] += 1
      end
      if p[8]==0
        u[52] += 1
      end
  elseif temp<cum_rate[3]
        u[8] += 1
  elseif temp<cum_rate[4]
      u[7] += 1
  elseif temp<cum_rate[5]
        u[6] += 1
  elseif temp<cum_rate[6]
        u[5] += 1
  elseif temp<cum_rate[7]
      if maximum([p[7],p[8]])==0
        u[1] -= 1
      end
        u[3] -=1
  elseif temp<cum_rate[8]
      if maximum([p[7],p[8]])==0
        u[2] -= 1
      end
        u[4] -=1
  elseif temp<cum_rate[9]
        u[7] -= 1
  elseif temp<cum_rate[10]
        u[8] -= 1
  elseif temp<cum_rate[11]
        u[5] -= 1
  elseif temp<cum_rate[12]
        u[6] -= 1
  elseif temp<cum_rate[13]
        u[49] -= 1
  elseif temp<cum_rate[14]
        u[50] -= 1
  elseif temp<cum_rate[15]
        u[51] -= 1
  elseif temp<cum_rate[16]
        u[52] -= 1
  elseif temp<cum_rate[17]
      if maximum([p[7],p[8]])>1
        u[10] += 1
        u[1] -= 1
      else
        u[1] -= 1
      end
      p[7]==1 ? u[49] += 1 : 0
      p[8]==1 ? u[50] += 1 : 0
  elseif temp<cum_rate[18]
      if maximum([p[7],p[8]])>2
        u[10] -= 1
        u[11] += 1
      else
        u[10] -= 1
      end
      p[7]==2 ? u[49] += 1 : 0
      p[8]==2 ? u[50] += 1 : 0
  elseif temp<cum_rate[19]
      if maximum([p[7],p[8]])>3
        u[11] -= 1
        u[12] += 1
      else
        u[11] -= 1
      end
      p[7]==3 ? u[49] += 1 : 0
      p[8]==3 ? u[50] += 1 : 0
  elseif temp<cum_rate[20]
      if maximum([p[7],p[8]])>4
        u[12] -= 1
        u[13] += 1
      else
        u[12] -= 1
      end
      p[7]==4 ? u[49] += 1 : 0
      p[8]==4 ? u[50] += 1 : 0
  elseif temp<cum_rate[21]
      if maximum([p[7],p[8]])>5
        u[13] -= 1
        u[14] += 1
      else
        u[13] -= 1
      end
      p[7]==5 ? u[49] += 1 : 0
      p[8]==5 ? u[50] += 1 : 0
  elseif temp<cum_rate[22]
      if maximum([p[7],p[8]])>6
        u[14] -= 1
        u[15] += 1
      else
        u[14] -= 1
      end
      p[7]==6 ? u[49] += 1 : 0
      p[8]==6 ? u[50] += 1 : 0
  elseif temp<cum_rate[23]
      if maximum([p[7],p[8]])>7
        u[15] -= 1
        u[16] += 1
      else
        u[15] -= 1
      end
      p[7]==7 ? u[49] += 1 : 0
      p[8]==7 ? u[50] += 1 : 0
  elseif temp<cum_rate[24]
      if maximum([p[7],p[8]])>8
        u[16] -= 1
        u[17] += 1
      else
        u[16] -= 1
      end
      p[7]==8 ? u[49] += 1 : 0
      p[8]==8 ? u[50] += 1 : 0
  elseif temp<cum_rate[25]
      if maximum([p[7],p[8]])>9
        u[17] -= 1
        u[18] += 1
      else
        u[17] -= 1
      end
      p[7]==9 ? u[49] += 1 : 0
      p[8]==9 ? u[50] += 1 : 0
  elseif temp<cum_rate[26]
      if maximum([p[7],p[8]])>10
        u[18] -= 1
        u[19] += 1
      else
        u[18] -= 1
      end
      p[7]==10 ? u[49] += 1 : 0
      p[8]==10 ? u[50] += 1 : 0
  elseif temp<cum_rate[27]
      if maximum([p[7],p[8]])>11
        u[19] -= 1
        u[20] += 1
      else
        u[19] -= 1
      end
      p[7]==11 ? u[49] += 1 : 0
      p[8]==11 ? u[50] += 1 : 0
  elseif temp<cum_rate[28]
      if maximum([p[7],p[8]])>12
        u[20] -= 1
        u[21] += 1
      else
        u[20] -= 1
      end
      p[7]==12 ? u[49] += 1 : 0
      p[8]==12 ? u[50] += 1 : 0
  elseif temp<cum_rate[29]
      if maximum([p[7],p[8]])>13
        u[21] -= 1
        u[22] += 1
      else
        u[21] -= 1
      end
      p[7]==13 ? u[49] += 1 : 0
      p[8]==13 ? u[50] += 1 : 0
  elseif temp<cum_rate[30]
      if maximum([p[7],p[8]])>14
        u[22] -= 1
        u[23] += 1
      else
        u[22] -= 1
      end
      p[7]==14 ? u[49] += 1 : 0
      p[8]==14 ? u[50] += 1 : 0
  elseif temp<cum_rate[31]
      if maximum([p[7],p[8]])>15
        u[23] -= 1
        u[24] += 1
      else
        u[23] -= 1
      end
      p[7]==15 ? u[49] += 1 : 0
      p[8]==15 ? u[50] += 1 : 0
  elseif temp<cum_rate[32]
      if maximum([p[7],p[8]])>16
        u[24] -= 1
        u[25] += 1
      else
        u[24] -= 1
      end
      p[7]==16 ? u[49] += 1 : 0
      p[8]==16 ? u[50] += 1 : 0
  elseif temp<cum_rate[33]
      if maximum([p[7],p[8]])>17
        u[25] -= 1
        u[26] += 1
      else
        u[25] -= 1
      end
      p[7]==17 ? u[49] += 1 : 0
      p[8]==17 ? u[50] += 1 : 0
  elseif temp<cum_rate[34]
      if maximum([p[7],p[8]])>18
        u[26] -= 1
        u[27] += 1
      else
        u[26] -= 1
      end
      p[7]==18 ? u[49] += 1 : 0
      p[8]==18 ? u[50] += 1 : 0
  elseif temp<cum_rate[35]
      if maximum([p[7],p[8]])>19
        u[27] -= 1
        u[28] += 1
      else
        u[27] -= 1
      end
      p[7]==19 ? u[49] += 1 : 0
      p[8]==19 ? u[50] += 1 : 0
  elseif temp<cum_rate[36]
      u[28] -= 1
      p[7]==20 ? u[49] += 1 : 0
      p[8]==20 ? u[50] += 1 : 0
  elseif temp<cum_rate[37]
      if maximum([p[7],p[8]])>1
        u[30] += 1
        u[2] -= 1
      else
        u[2] -=1
      end
      p[7]==1 ? u[51] += 1 : 0
      p[8]==1 ? u[52] += 1 : 0
  elseif temp<cum_rate[38]
      if maximum([p[7],p[8]])>2
        u[30] -= 1
        u[31] += 1
      else
        u[30] -= 1
      end
      p[7]==2 ? u[51] += 1 : 0
      p[8]==2 ? u[52] += 1 : 0
  elseif temp<cum_rate[39]
      if maximum([p[7],p[8]])>3
        u[31] -= 1
        u[32] += 1
      else
        u[31] -= 1
      end
      p[7]==3 ? u[51] += 1 : 0
      p[8]==3 ? u[52] += 1 : 0
  elseif temp<cum_rate[40]
      if maximum([p[7],p[8]])>4
        u[32] -= 1
        u[33] += 1
      else
        u[32] -= 1
      end
      p[7]==4 ? u[51] += 1 : 0
      p[8]==4 ? u[52] += 1 : 0
  elseif temp<cum_rate[41]
      if maximum([p[7],p[8]])>5
        u[33] -= 1
        u[34] += 1
      else
        u[33] -= 1
      end
      p[7]==5 ? u[51] += 1 : 0
      p[8]==5 ? u[52] += 1 : 0
  elseif temp<cum_rate[42]
      if maximum([p[7],p[8]])>6
        u[34] -= 1
        u[35] += 1
      else
        u[34] -= 1
      end
      p[7]==6 ? u[51] += 1 : 0
      p[8]==6 ? u[52] += 1 : 0
  elseif temp<cum_rate[43]
      if maximum([p[7],p[8]])>7
        u[35] -= 1
        u[36] += 1
      else
        u[35] -= 1
      end
      p[7]==7 ? u[51] += 1 : 0
      p[8]==7 ? u[52] += 1 : 0
  elseif temp<cum_rate[44]
      if maximum([p[7],p[8]])>8
        u[36] -= 1
        u[37] += 1
      else
        u[36] -= 1
      end
      p[7]==8 ? u[51] += 1 : 0
      p[8]==8 ? u[52] += 1 : 0
  elseif temp<cum_rate[45]
      if maximum([p[7],p[8]])>9
        u[37] -= 1
        u[38] += 1
      else
        u[37] -= 1
      end
      p[7]==9 ? u[51] += 1 : 0
      p[8]==9 ? u[52] += 1 : 0
  elseif temp<cum_rate[46]
      if maximum([p[7],p[8]])>10
        u[38] -= 1
        u[39] += 1
      else
        u[38] -= 1
      end
      p[7]==10 ? u[51] += 1 : 0
      p[8]==10 ? u[52] += 1 : 0
  elseif temp<cum_rate[47]
      if maximum([p[7],p[8]])>11
        u[39] -= 1
        u[40] += 1
      else
        u[39] -= 1
      end
      p[7]==11 ? u[51] += 1 : 0
      p[8]==11 ? u[52] += 1 : 0
  elseif temp<cum_rate[48]
      if maximum([p[7],p[8]])>12
        u[40] -= 1
        u[41] += 1
      else
        u[40] -= 1
      end
      p[7]==12 ? u[51] += 1 : 0
      p[8]==12 ? u[52] += 1 : 0
  elseif temp<cum_rate[49]
      if maximum([p[7],p[8]])>13
        u[41] -= 1
        u[42] += 1
      else
        u[41] -= 1
      end
      p[7]==13 ? u[51] += 1 : 0
      p[8]==13 ? u[52] += 1 : 0
  elseif temp<cum_rate[50]
      if maximum([p[7],p[8]])>14
        u[42] -= 1
        u[43] += 1
      else
        u[42] -= 1
      end
      p[7]==14 ? u[51] += 1 : 0
      p[8]==14 ? u[52] += 1 : 0
  elseif temp<cum_rate[51]
      if maximum([p[7],p[8]])>15
        u[43] -= 1
        u[44] += 1
      else
        u[43] -= 1
      end
      p[7]==15 ? u[51] += 1 : 0
      p[8]==15 ? u[52] += 1 : 0
  elseif temp<cum_rate[52]
      if maximum([p[7],p[8]])>16
        u[44] -= 1
        u[45] += 1
      else
        u[44] -= 1
      end
      p[7]==16 ? u[51] += 1 : 0
      p[8]==16 ? u[52] += 1 : 0
  elseif temp<cum_rate[53]
      if maximum([p[7],p[8]])>17
        u[45] -= 1
        u[46] += 1
      else
        u[45] -= 1
      end
      p[7]==17 ? u[51] += 1 : 0
      p[8]==17 ? u[52] += 1 : 0
  elseif temp<cum_rate[54]
      if maximum([p[7],p[8]])>18
        u[46] -= 1
        u[47] += 1
      else
        u[46] -= 1
      end
      p[7]==18 ? u[51] += 1 : 0
      p[8]==18 ? u[52] += 1 : 0
  elseif temp<cum_rate[55]
      if maximum([p[7],p[8]])>19
        u[47] -= 1
        u[48] += 1
      else
        u[47] -= 1
      end
      p[7]==19 ? u[51] += 1 : 0
      p[8]==19 ? u[52] += 1 : 0
  elseif temp<cum_rate[56]
      u[48] -= 1
      p[7]==20 ? u[51] += 1 : 0
      p[8]==20 ? u[52] += 1 : 0
  end
    t = t+t_step
    push!(out,copy(u))
    push!(t_out,copy(t))
  end
  return out, t_out
end

### Aneuploid cells
function sim_delay_aneu(p::Array{Float64,1},tspan::Tuple{Float64,Float64},u0::Array{Int64,1})
  t::Float64=tspan[1]
  i::Int64=0
  out::Array{Array{Int64,1}} = []
  t_out::Array{Float64,1} = []
  push!(out, u0)
  push!(t_out,t)
  u::Array{Int64,1}=copy(u0)
  #############################################################
  rate::Array{Float64,1}=zeros(56)
  #more than 1X?
  p1X::Int64 = 0
  if p[26] == 3
      rate = [rate; zeros(28)]
  elseif p[26] ==4
      rate = [rate; zeros(56)]
  end
  p[24:25] = 0
  p[27:32] = 0

  if p[26]>=1
      p[29] = 1
      if p[26] >=2
          p1X = 1
          p[30] = 1

          if p[26]>=3
              p[31] = 1
              if p[26]>=4
                  p[32] = 1
              end
          end
      end
  end
  while (t<tspan[2])
    i+=1
    if p[26] <= 2
        # Xist transcription
        rate[1] = p[21]*(p[24]+p[29]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[7]^p[13]/(u[7]^p[13]+(p[22]*p[14])^p[13])))
        rate[2] = p1X*p[21]*(p[25]+p[30]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[8]^p[13]/(u[8]^p[13]+(p[22]*p[14])^p[13])))
        # cXR production
        rate[3] = p1X*p[22]*(1-u[51]^p[5]/(u[51]^p[5]+(p[21]*p[6])^p[5]))
        rate[4] = p[22]*(1-u[49]^p[5]/(u[49]^p[5]+(p[21]*p[6])^p[5]))
        # tXA production
        rate[5] = p1X*p[23]*(1-u[52]^p[3]/(u[52]^p[3]+(p[21]*p[4])^p[3]))
        rate[6] = p[23]*(1-u[50]^p[3]/(u[50]^p[3]+(p[21]*p[4])^p[3]))
        # degradation
        rate[7] = 0.1733*u[3]
        rate[8] = 0.1733*u[4]
        rate[9] = 1*u[7]
        rate[10] = 1*u[8]
        rate[11] = 1*u[5]
        rate[12] = 1*u[6]
        rate[13] = 1*u[49]
        rate[14] = 1*u[50]
        rate[15] = 1*u[51]
        rate[16] = 1*u[52]
        ### silencing intermediates on chr1
        rate[17] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[1],0)
        rate[18] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[10],0)
        rate[19] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[11],0)
        rate[20] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[12],0)
        rate[21] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[13],0)
        rate[22] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[14],0)
        rate[23] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[15],0)
        rate[24] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[16],0)
        rate[25] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[17],0)
        rate[26] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[18],0)
        rate[27] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[19],0)
        rate[28] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[20],0)
        rate[29] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[21],0)
        rate[30] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[22],0)
        rate[31] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[23],0)
        rate[32] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[24],0)
        rate[33] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[25],0)
        rate[34] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[26],0)
        rate[35] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[27],0)
        rate[36] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[28],0)
        ### silencing intermediates chr 2
        rate[37] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[2],0)
        rate[38] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[30],0)
        rate[39] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[31],0)
        rate[40] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[32],0)
        rate[41] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[33],0)
        rate[42] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[34],0)
        rate[43] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[35],0)
        rate[44] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[36],0)
        rate[45] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[37],0)
        rate[46] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[38],0)
        rate[47] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[39],0)
        rate[48] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[40],0)
        rate[49] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[41],0)
        rate[50] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[42],0)
        rate[51] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[43],0)
        rate[52] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[44],0)
        rate[53] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[45],0)
        rate[54] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[46],0)
        rate[55] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[47],0)
        rate[56] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[48],0)
    elseif p[26] == 3
        # Xist transcription
        rate[1] = p[21]*(p[24]+p[29]*(0.5*(u[5]+u[6]+u[54]))^p[11]/((0.5*(u[5]+u[6]+u[54]))^p[11]+(p[23]*p[12])^p[11])*(1-u[7]^p[13]/(u[7]^p[13]+(p[22]*p[14])^p[13])))
        rate[2] = p[21]*(p[25]+p[30]*(0.5*(u[5]+u[6]+u[54]))^p[11]/((0.5*(u[5]+u[6]+u[54]))^p[11]+(p[23]*p[12])^p[11])*(1-u[8]^p[13]/(u[8]^p[13]+(p[22]*p[14])^p[13])))
        rate[57] = p[21]*(p[27]+p[31]*(0.5*(u[5]+u[6]+u[54]))^p[11]/((0.5*(u[5]+u[6]+u[54]))^p[11]+(p[23]*p[12])^p[11])*(1-u[55]^p[13]/(u[55]^p[13]+(p[22]*p[14])^p[13])))
        # cXR production (scaling factor * (1-silencing))
        rate[3] = p[22]*(1-u[51]^p[5]/(u[51]^p[5]+(p[21]*p[6])^p[5]))
        rate[4] = p[22]*(1-u[49]^p[5]/(u[49]^p[5]+(p[21]*p[6])^p[5]))
        rate[58] = p[22]*(1-u[97]^p[5]/(u[97]^p[5]+(p[21]*p[6])^p[5]))
        # tXA production (scaling factor * (1-silencing))
        rate[5] = p[23]*(1-u[52]^p[3]/(u[52]^p[3]+(p[21]*p[4])^p[3]))
        rate[6] = p[23]*(1-u[50]^p[3]/(u[50]^p[3]+(p[21]*p[4])^p[3]))
        rate[59] = p[23]*(1-u[98]^p[3]/(u[98]^p[3]+(p[21]*p[4])^p[3]))
        # degradation
        rate[7] = 0.1733*u[3]
        rate[8] = 0.1733*u[4]
        rate[60] = 0.1733*u[101]
        rate[9] = 1*u[7]
        rate[10] = 1*u[8]
        rate[61] = 1*u[55]
        rate[11] = 1*u[5]
        rate[12] = 1*u[6]
        rate[62] = 1*u[54]
        rate[13] = 1*u[49]
        rate[14] = 1*u[50]
        rate[63] = 1*u[97]
        rate[15] = 1*u[51]
        rate[16] = 1*u[52]
        rate[64] = 1*u[98]
        ### silencing intermediates on chr1
        rate[17] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[1],0)
        rate[18] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[10],0)
        rate[19] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[11],0)
        rate[20] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[12],0)
        rate[21] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[13],0)
        rate[22] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[14],0)
        rate[23] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[15],0)
        rate[24] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[16],0)
        rate[25] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[17],0)
        rate[26] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[18],0)
        rate[27] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[19],0)
        rate[28] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[20],0)
        rate[29] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[21],0)
        rate[30] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[22],0)
        rate[31] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[23],0)
        rate[32] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[24],0)
        rate[33] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[25],0)
        rate[34] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[26],0)
        rate[35] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[27],0)
        rate[36] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[28],0)
        ### silencing intermediates chr 2
        rate[37] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[2],0)
        rate[38] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[30],0)
        rate[39] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[31],0)
        rate[40] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[32],0)
        rate[41] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[33],0)
        rate[42] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[34],0)
        rate[43] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[35],0)
        rate[44] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[36],0)
        rate[45] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[37],0)
        rate[46] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[38],0)
        rate[47] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[39],0)
        rate[48] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[40],0)
        rate[49] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[41],0)
        rate[50] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[42],0)
        rate[51] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[43],0)
        rate[52] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[44],0)
        rate[53] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[45],0)
        rate[54] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[46],0)
        rate[55] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[47],0)
        rate[56] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[48],0)
        ### silencing intermediates chr 3
        rate[65] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[53],0)
        rate[66] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[59],0)
        rate[67] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[60],0)
        rate[68] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[61],0)
        rate[69] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[62],0)
        rate[70] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[63],0)
        rate[71] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[64],0)
        rate[72] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[65],0)
        rate[73] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[66],0)
        rate[74] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[67],0)
        rate[75] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[68],0)
        rate[76] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[69],0)
        rate[77] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[70],0)
        rate[78] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[71],0)
        rate[79] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[72],0)
        rate[80] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[73],0)
        rate[81] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[74],0)
        rate[82] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[75],0)
        rate[83] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[76],0)
        rate[84] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[77],0)
    elseif p[26] == 4
        # Xist transcription
        rate[1] = p[21]*(p[24]+p[29]*(0.5*(u[5]+u[6]+u[54]+u[57]))^p[11]/((0.5*(u[5]+u[6]+u[54]+u[57]))^p[11]+(p[23]*p[12])^p[11])*(1-u[7]^p[13]/(u[7]^p[13]+(p[22]*p[14])^p[13])))
        rate[2] = p[21]*(p[25]+p[30]*(0.5*(u[5]+u[6]+u[54]+u[57]))^p[11]/((0.5*(u[5]+u[6]+u[54]+u[57]))^p[11]+(p[23]*p[12])^p[11])*(1-u[8]^p[13]/(u[8]^p[13]+(p[22]*p[14])^p[13])))
        rate[57] = p[21]*(p[27]+p[31]*(0.5*(u[5]+u[6]+u[54]+u[57]))^p[11]/((0.5*(u[5]+u[6]+u[54]+u[57]))^p[11]+(p[23]*p[12])^p[11])*(1-u[55]^p[13]/(u[55]^p[13]+(p[22]*p[14])^p[13])))
        rate[85] = p[21]*(p[28]+p[32]*(0.5*(u[5]+u[6]+u[54]+u[57]))^p[11]/((0.5*(u[5]+u[6]+u[54]+u[57]))^p[11]+(p[23]*p[12])^p[11])*(1-u[58]^p[13]/(u[58]^p[13]+(p[22]*p[14])^p[13])))
        # cXR production (scaling factor * (1-silencing))
        rate[3] = p[22]*(1-u[51]^p[5]/(u[51]^p[5]+(p[21]*p[6])^p[5]))
        rate[4] = p[22]*(1-u[49]^p[5]/(u[49]^p[5]+(p[21]*p[6])^p[5]))
        rate[58] = p[22]*(1-u[97]^p[5]/(u[97]^p[5]+(p[21]*p[6])^p[5]))
        rate[86] = p[22]*(1-u[99]^p[5]/(u[99]^p[5]+(p[21]*p[6])^p[5]))
        # tXA production (scaling factor * (1-silencing))
        rate[5] = p[23]*(1-u[52]^p[3]/(u[52]^p[3]+(p[21]*p[4])^p[3]))
        rate[6] = p[23]*(1-u[50]^p[3]/(u[50]^p[3]+(p[21]*p[4])^p[3]))
        rate[59] = p[23]*(1-u[98]^p[3]/(u[98]^p[3]+(p[21]*p[4])^p[3]))
        rate[87] = p[23]*(1-u[100]^p[3]/(u[100]^p[3]+(p[21]*p[4])^p[3]))
        # degradation
        rate[7] = 0.1733*u[3]
        rate[8] = 0.1733*u[4]
        rate[60] = 0.1733*u[101]
        rate[88] = 0.1733*u[102]
        rate[9] = 1*u[7]
        rate[10] = 1*u[8]
        rate[61] = 1*u[55]
        rate[89] = 1*u[58]
        rate[11] = 1*u[5]
        rate[12] = 1*u[6]
        rate[62] = 1*u[54]
        rate[90] = 1*u[57]
        rate[13] = 1*u[49]
        rate[14] = 1*u[50]
        rate[63] = 1*u[97]
        rate[91] = 1*u[99]
        rate[15] = 1*u[51]
        rate[16] = 1*u[52]
        rate[64] = 1*u[98]
        rate[92] = 1*u[100]
        ### silencing intermediates on chr1
        rate[17] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[1],0)
        rate[18] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[10],0)
        rate[19] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[11],0)
        rate[20] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[12],0)
        rate[21] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[13],0)
        rate[22] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[14],0)
        rate[23] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[15],0)
        rate[24] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[16],0)
        rate[25] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[17],0)
        rate[26] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[18],0)
        rate[27] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[19],0)
        rate[28] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[20],0)
        rate[29] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[21],0)
        rate[30] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[22],0)
        rate[31] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[23],0)
        rate[32] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[24],0)
        rate[33] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[25],0)
        rate[34] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[26],0)
        rate[35] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[27],0)
        rate[36] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[28],0)
        ### silencing intermediates chr 2
        rate[37] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[2],0)
        rate[38] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[30],0)
        rate[39] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[31],0)
        rate[40] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[32],0)
        rate[41] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[33],0)
        rate[42] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[34],0)
        rate[43] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[35],0)
        rate[44] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[36],0)
        rate[45] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[37],0)
        rate[46] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[38],0)
        rate[47] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[39],0)
        rate[48] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[40],0)
        rate[49] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[41],0)
        rate[50] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[42],0)
        rate[51] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[43],0)
        rate[52] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[44],0)
        rate[53] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[45],0)
        rate[54] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[46],0)
        rate[55] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[47],0)
        rate[56] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[48],0)
        ### silencing intermediates chr 3
        rate[65] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[53],0)
        rate[66] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[59],0)
        rate[67] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[60],0)
        rate[68] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[61],0)
        rate[69] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[62],0)
        rate[70] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[63],0)
        rate[71] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[64],0)
        rate[72] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[65],0)
        rate[73] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[66],0)
        rate[74] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[67],0)
        rate[75] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[68],0)
        rate[76] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[69],0)
        rate[77] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[70],0)
        rate[78] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[71],0)
        rate[79] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[72],0)
        rate[80] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[73],0)
        rate[81] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[74],0)
        rate[82] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[75],0)
        rate[83] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[76],0)
        rate[84] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[77],0)
        ### silencing intermediates chr 4
        rate[93] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[56],0)
        rate[94] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[78],0)
        rate[95] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[79],0)
        rate[96] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[80],0)
        rate[97] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[81],0)
        rate[98] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[82],0)
        rate[99] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[83],0)
        rate[100] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[84],0)
        rate[101] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[85],0)
        rate[102] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[86],0)
        rate[103] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[87],0)
        rate[104] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[88],0)
        rate[105] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[89],0)
        rate[106] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[90],0)
        rate[107] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[91],0)
        rate[108] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[92],0)
        rate[109] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[93],0)
        rate[110] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[94],0)
        rate[111] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[95],0)
        rate[112] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[96],0)
    end
	
	# Gillespie
    r = rand(2)
    # determine time point of next reaction
    t_step::Float64 = 1/sum(rate)*log(1/(1-r[1]))
    # decide which reaction will occur
    temp::Float64 = r[2]*sum(rate)
    cum_rate = cumsum(rate)
    # Execute reaction
    if temp<cum_rate[1]
      u[1] += 1
      u[3] += 1
        if p[7]==0
        u[49] += 1
      end
      if p[8]==0
        u[50] += 1
      end
  elseif temp<cum_rate[2]
      u[2] += 1
      u[4] += 1
      if p[7]==0
        u[51] += 1
      end
      if p[8]==0
        u[52] += 1
      end
  elseif temp<cum_rate[3]
        u[8] += 1
  elseif temp<cum_rate[4]
      u[7] += 1
  elseif temp<cum_rate[5]
        u[6] += 1
  elseif temp<cum_rate[6]
        u[5] += 1
  elseif temp<cum_rate[7]
      if maximum([p[7],p[8]])==0
        u[1] -= 1
      end
        u[3] -=1
  elseif temp<cum_rate[8]
      if maximum([p[7],p[8]])==0
        u[2] -= 1
      end
        u[4] -=1
  elseif temp<cum_rate[9]
        u[7] -= 1
  elseif temp<cum_rate[10]
        u[8] -= 1
  elseif temp<cum_rate[11]
        u[5] -= 1
  elseif temp<cum_rate[12]
        u[6] -= 1
  elseif temp<cum_rate[13]
        u[49] -= 1
  elseif temp<cum_rate[14]
        u[50] -= 1
  elseif temp<cum_rate[15]
        u[51] -= 1
  elseif temp<cum_rate[16]
        u[52] -= 1
  elseif temp<cum_rate[17]
      if maximum([p[7],p[8]])>1
        u[10] += 1
        u[1] -= 1
      else
        u[1] -= 1
      end
      p[7]==1 ? u[49] += 1 : 0
      p[8]==1 ? u[50] += 1 : 0
  elseif temp<cum_rate[18]
      if maximum([p[7],p[8]])>2
        u[10] -= 1
        u[11] += 1
      else
        u[10] -= 1
      end
      p[7]==2 ? u[49] += 1 : 0
      p[8]==2 ? u[50] += 1 : 0
  elseif temp<cum_rate[19]
      if maximum([p[7],p[8]])>3
        u[11] -= 1
        u[12] += 1
      else
        u[11] -= 1
      end
      p[7]==3 ? u[49] += 1 : 0
      p[8]==3 ? u[50] += 1 : 0
  elseif temp<cum_rate[20]
      if maximum([p[7],p[8]])>4
        u[12] -= 1
        u[13] += 1
      else
        u[12] -= 1
      end
      p[7]==4 ? u[49] += 1 : 0
      p[8]==4 ? u[50] += 1 : 0
  elseif temp<cum_rate[21]
      if maximum([p[7],p[8]])>5
        u[13] -= 1
        u[14] += 1
      else
        u[13] -= 1
      end
      p[7]==5 ? u[49] += 1 : 0
      p[8]==5 ? u[50] += 1 : 0
  elseif temp<cum_rate[22]
      if maximum([p[7],p[8]])>6
        u[14] -= 1
        u[15] += 1
      else
        u[14] -= 1
      end
      p[7]==6 ? u[49] += 1 : 0
      p[8]==6 ? u[50] += 1 : 0
  elseif temp<cum_rate[23]
      if maximum([p[7],p[8]])>7
        u[15] -= 1
        u[16] += 1
      else
        u[15] -= 1
      end
      p[7]==7 ? u[49] += 1 : 0
      p[8]==7 ? u[50] += 1 : 0
  elseif temp<cum_rate[24]
      if maximum([p[7],p[8]])>8
        u[16] -= 1
        u[17] += 1
      else
        u[16] -= 1
      end
      p[7]==8 ? u[49] += 1 : 0
      p[8]==8 ? u[50] += 1 : 0
  elseif temp<cum_rate[25]
      if maximum([p[7],p[8]])>9
        u[17] -= 1
        u[18] += 1
      else
        u[17] -= 1
      end
      p[7]==9 ? u[49] += 1 : 0
      p[8]==9 ? u[50] += 1 : 0
  elseif temp<cum_rate[26]
      if maximum([p[7],p[8]])>10
        u[18] -= 1
        u[19] += 1
      else
        u[18] -= 1
      end
      p[7]==10 ? u[49] += 1 : 0
      p[8]==10 ? u[50] += 1 : 0
  elseif temp<cum_rate[27]
      if maximum([p[7],p[8]])>11
        u[19] -= 1
        u[20] += 1
      else
        u[19] -= 1
      end
      p[7]==11 ? u[49] += 1 : 0
      p[8]==11 ? u[50] += 1 : 0
  elseif temp<cum_rate[28]
      if maximum([p[7],p[8]])>12
        u[20] -= 1
        u[21] += 1
      else
        u[20] -= 1
      end
      p[7]==12 ? u[49] += 1 : 0
      p[8]==12 ? u[50] += 1 : 0
  elseif temp<cum_rate[29]
      if maximum([p[7],p[8]])>13
        u[21] -= 1
        u[22] += 1
      else
        u[21] -= 1
      end
      p[7]==13 ? u[49] += 1 : 0
      p[8]==13 ? u[50] += 1 : 0
  elseif temp<cum_rate[30]
      if maximum([p[7],p[8]])>14
        u[22] -= 1
        u[23] += 1
      else
        u[22] -= 1
      end
      p[7]==14 ? u[49] += 1 : 0
      p[8]==14 ? u[50] += 1 : 0
  elseif temp<cum_rate[31]
      if maximum([p[7],p[8]])>15
        u[23] -= 1
        u[24] += 1
      else
        u[23] -= 1
      end
      p[7]==15 ? u[49] += 1 : 0
      p[8]==15 ? u[50] += 1 : 0
  elseif temp<cum_rate[32]
      if maximum([p[7],p[8]])>16
        u[24] -= 1
        u[25] += 1
      else
        u[24] -= 1
      end
      p[7]==16 ? u[49] += 1 : 0
      p[8]==16 ? u[50] += 1 : 0
  elseif temp<cum_rate[33]
      if maximum([p[7],p[8]])>17
        u[25] -= 1
        u[26] += 1
      else
        u[25] -= 1
      end
      p[7]==17 ? u[49] += 1 : 0
      p[8]==17 ? u[50] += 1 : 0
  elseif temp<cum_rate[34]
      if maximum([p[7],p[8]])>18
        u[26] -= 1
        u[27] += 1
      else
        u[26] -= 1
      end
      p[7]==18 ? u[49] += 1 : 0
      p[8]==18 ? u[50] += 1 : 0
  elseif temp<cum_rate[35]
      if maximum([p[7],p[8]])>19
        u[27] -= 1
        u[28] += 1
      else
        u[27] -= 1
      end
      p[7]==19 ? u[49] += 1 : 0
      p[8]==19 ? u[50] += 1 : 0
  elseif temp<cum_rate[36]
      u[28] -= 1
      p[7]==20 ? u[49] += 1 : 0
      p[8]==20 ? u[50] += 1 : 0
  elseif temp<cum_rate[37]
      if maximum([p[7],p[8]])>1
        u[30] += 1
        u[2] -= 1
      else
        u[2] -=1
      end
      p[7]==1 ? u[51] += 1 : 0
      p[8]==1 ? u[52] += 1 : 0
  elseif temp<cum_rate[38]
      if maximum([p[7],p[8]])>2
        u[30] -= 1
        u[31] += 1
      else
        u[30] -= 1
      end
      p[7]==2 ? u[51] += 1 : 0
      p[8]==2 ? u[52] += 1 : 0
  elseif temp<cum_rate[39]
      if maximum([p[7],p[8]])>3
        u[31] -= 1
        u[32] += 1
      else
        u[31] -= 1
      end
      p[7]==3 ? u[51] += 1 : 0
      p[8]==3 ? u[52] += 1 : 0
  elseif temp<cum_rate[40]
      if maximum([p[7],p[8]])>4
        u[32] -= 1
        u[33] += 1
      else
        u[32] -= 1
      end
      p[7]==4 ? u[51] += 1 : 0
      p[8]==4 ? u[52] += 1 : 0
  elseif temp<cum_rate[41]
      if maximum([p[7],p[8]])>5
        u[33] -= 1
        u[34] += 1
      else
        u[33] -= 1
      end
      p[7]==5 ? u[51] += 1 : 0
      p[8]==5 ? u[52] += 1 : 0
  elseif temp<cum_rate[42]
      if maximum([p[7],p[8]])>6
        u[34] -= 1
        u[35] += 1
      else
        u[34] -= 1
      end
      p[7]==6 ? u[51] += 1 : 0
      p[8]==6 ? u[52] += 1 : 0
  elseif temp<cum_rate[43]
      if maximum([p[7],p[8]])>7
        u[35] -= 1
        u[36] += 1
      else
        u[35] -= 1
      end
      p[7]==7 ? u[51] += 1 : 0
      p[8]==7 ? u[52] += 1 : 0
  elseif temp<cum_rate[44]
      if maximum([p[7],p[8]])>8
        u[36] -= 1
        u[37] += 1
      else
        u[36] -= 1
      end
      p[7]==8 ? u[51] += 1 : 0
      p[8]==8 ? u[52] += 1 : 0
  elseif temp<cum_rate[45]
      if maximum([p[7],p[8]])>9
        u[37] -= 1
        u[38] += 1
      else
        u[37] -= 1
      end
      p[7]==9 ? u[51] += 1 : 0
      p[8]==9 ? u[52] += 1 : 0
  elseif temp<cum_rate[46]
      if maximum([p[7],p[8]])>10
        u[38] -= 1
        u[39] += 1
      else
        u[38] -= 1
      end
      p[7]==10 ? u[51] += 1 : 0
      p[8]==10 ? u[52] += 1 : 0
  elseif temp<cum_rate[47]
      if maximum([p[7],p[8]])>11
        u[39] -= 1
        u[40] += 1
      else
        u[39] -= 1
      end
      p[7]==11 ? u[51] += 1 : 0
      p[8]==11 ? u[52] += 1 : 0
  elseif temp<cum_rate[48]
      if maximum([p[7],p[8]])>12
        u[40] -= 1
        u[41] += 1
      else
        u[40] -= 1
      end
      p[7]==12 ? u[51] += 1 : 0
      p[8]==12 ? u[52] += 1 : 0
  elseif temp<cum_rate[49]
      if maximum([p[7],p[8]])>13
        u[41] -= 1
        u[42] += 1
      else
        u[41] -= 1
      end
      p[7]==13 ? u[51] += 1 : 0
      p[8]==13 ? u[52] += 1 : 0
  elseif temp<cum_rate[50]
      if maximum([p[7],p[8]])>14
        u[42] -= 1
        u[43] += 1
      else
        u[42] -= 1
      end
      p[7]==14 ? u[51] += 1 : 0
      p[8]==14 ? u[52] += 1 : 0
  elseif temp<cum_rate[51]
      if maximum([p[7],p[8]])>15
        u[43] -= 1
        u[44] += 1
      else
        u[43] -= 1
      end
      p[7]==15 ? u[51] += 1 : 0
      p[8]==15 ? u[52] += 1 : 0
  elseif temp<cum_rate[52]
      if maximum([p[7],p[8]])>16
        u[44] -= 1
        u[45] += 1
      else
        u[44] -= 1
      end
      p[7]==16 ? u[51] += 1 : 0
      p[8]==16 ? u[52] += 1 : 0
  elseif temp<cum_rate[53]
      if maximum([p[7],p[8]])>17
        u[45] -= 1
        u[46] += 1
      else
        u[45] -= 1
      end
      p[7]==17 ? u[51] += 1 : 0
      p[8]==17 ? u[52] += 1 : 0
  elseif temp<cum_rate[54]
      if maximum([p[7],p[8]])>18
        u[46] -= 1
        u[47] += 1
      else
        u[46] -= 1
      end
      p[7]==18 ? u[51] += 1 : 0
      p[8]==18 ? u[52] += 1 : 0
  elseif temp<cum_rate[55]
      if maximum([p[7],p[8]])>19
        u[47] -= 1
        u[48] += 1
      else
        u[47] -= 1
      end
      p[7]==19 ? u[51] += 1 : 0
      p[8]==19 ? u[52] += 1 : 0
  elseif temp<cum_rate[56]
      u[48] -= 1
      p[7]==20 ? u[51] += 1 : 0
      p[8]==20 ? u[52] += 1 : 0
  elseif temp<cum_rate[57]
      u[53] += 1
      u[101] += 1
      if p[7]==0
        u[97] += 1
      end
      if p[8]==0
        u[98] += 1
      end
  elseif temp<cum_rate[58]
      u[55] += 1
  elseif temp<cum_rate[59]
      u[54] += 1
  elseif temp<cum_rate[60]
      if maximum([p[7],p[8]])==0
        u[53] -= 1
      end
        u[101] -=1
  elseif temp<cum_rate[61]
      u[55] -= 1
  elseif temp<cum_rate[62]
      u[54] -= 1
  elseif temp<cum_rate[63]
      u[97] -= 1
  elseif temp<cum_rate[64]
      u[98] -= 1
  elseif temp<cum_rate[65]
      if maximum([p[7],p[8]])>1
        u[59] += 1
        u[53] -= 1
      else
        u[53] -=1
      end
      p[7]==1 ? u[97] += 1 : 0
      p[8]==1 ? u[98] += 1 : 0
  elseif temp<cum_rate[66]
      if maximum([p[7],p[8]])>2
        u[59] -= 1
        u[60] += 1
      else
        u[59] -= 1
      end
      p[7]==2 ? u[97] += 1 : 0
      p[8]==2 ? u[98] += 1 : 0
  elseif temp<cum_rate[67]
      if maximum([p[7],p[8]])>3
        u[60] -= 1
        u[61] += 1
      else
        u[60] -= 1
      end
      p[7]==3 ? u[97] += 1 : 0
      p[8]==3 ? u[98] += 1 : 0
  elseif temp<cum_rate[68]
      if maximum([p[7],p[8]])>4
        u[61] -= 1
        u[62] += 1
      else
        u[61] -= 1
      end
      p[7]==4 ? u[97] += 1 : 0
      p[8]==4 ? u[98] += 1 : 0
  elseif temp<cum_rate[69]
      if maximum([p[7],p[8]])>5
        u[62] -= 1
        u[63] += 1
      else
        u[62] -= 1
      end
      p[7]==5 ? u[97] += 1 : 0
      p[8]==5 ? u[98] += 1 : 0
  elseif temp<cum_rate[70]
      if maximum([p[7],p[8]])>6
        u[63] -= 1
        u[64] += 1
      else
        u[63] -= 1
      end
      p[7]==6 ? u[97] += 1 : 0
      p[8]==6 ? u[98] += 1 : 0
  elseif temp<cum_rate[71]
      if maximum([p[7],p[8]])>7
        u[64] -= 1
        u[65] += 1
      else
        u[64] -= 1
      end
      p[7]==7 ? u[97] += 1 : 0
      p[8]==7 ? u[98] += 1 : 0
  elseif temp<cum_rate[72]
      if maximum([p[7],p[8]])>8
        u[65] -= 1
        u[66] += 1
      else
        u[65] -= 1
      end
      p[7]==8 ? u[97] += 1 : 0
      p[8]==8 ? u[98] += 1 : 0
  elseif temp<cum_rate[73]
      if maximum([p[7],p[8]])>9
        u[66] -= 1
        u[67] += 1
      else
        u[66] -= 1
      end
      p[7]==9 ? u[97] += 1 : 0
      p[8]==9 ? u[98] += 1 : 0
  elseif temp<cum_rate[74]
      if maximum([p[7],p[8]])>10
        u[67] -= 1
        u[68] += 1
      else
        u[67] -= 1
      end
      p[7]==10 ? u[97] += 1 : 0
      p[8]==10 ? u[98] += 1 : 0
  elseif temp<cum_rate[75]
      if maximum([p[7],p[8]])>11
        u[68] -= 1
        u[69] += 1
      else
        u[68] -= 1
      end
      p[7]==11 ? u[97] += 1 : 0
      p[8]==11 ? u[98] += 1 : 0
  elseif temp<cum_rate[76]
      if maximum([p[7],p[8]])>12
        u[69] -= 1
        u[70] += 1
      else
        u[69] -= 1
      end
      p[7]==12 ? u[97] += 1 : 0
      p[8]==12 ? u[98] += 1 : 0
  elseif temp<cum_rate[77]
      if maximum([p[7],p[8]])>13
        u[70] -= 1
        u[71] += 1
      else
        u[70] -= 1
      end
      p[7]==13 ? u[97] += 1 : 0
      p[8]==13 ? u[98] += 1 : 0
  elseif temp<cum_rate[78]
      if maximum([p[7],p[8]])>14
        u[71] -= 1
        u[72] += 1
      else
        u[71] -= 1
      end
      p[7]==14 ? u[97] += 1 : 0
      p[8]==14 ? u[98] += 1 : 0
  elseif temp<cum_rate[79]
      if maximum([p[7],p[8]])>15
        u[72] -= 1
        u[73] += 1
      else
        u[72] -= 1
      end
      p[7]==15 ? u[97] += 1 : 0
      p[8]==15 ? u[98] += 1 : 0
  elseif temp<cum_rate[80]
      if maximum([p[7],p[8]])>16
        u[73] -= 1
        u[74] += 1
      else
        u[73] -= 1
      end
      p[7]==16 ? u[97] += 1 : 0
      p[8]==16 ? u[98] += 1 : 0
  elseif temp<cum_rate[81]
      if maximum([p[7],p[8]])>17
        u[74] -= 1
        u[75] += 1
      else
        u[74] -= 1
      end
      p[7]==17 ? u[97] += 1 : 0
      p[8]==17 ? u[98] += 1 : 0
  elseif temp<cum_rate[82]
      if maximum([p[7],p[8]])>18
        u[75] -= 1
        u[76] += 1
      else
        u[75] -= 1
      end
      p[7]==18 ? u[97] += 1 : 0
      p[8]==18 ? u[98] += 1 : 0
  elseif temp<cum_rate[83]
      if maximum([p[7],p[8]])>19
        u[76] -= 1
        u[77] += 1
      else
        u[76] -= 1
      end
      p[7]==19 ? u[97] += 1 : 0
      p[8]==19 ? u[98] += 1 : 0
  elseif temp<cum_rate[84]
      u[77] -= 1
      p[7]==20 ? u[97] += 1 : 0
      p[8]==20 ? u[98] += 1 : 0
  elseif temp<cum_rate[85]
      u[56] += 1
      u[102] += 1
      if p[7]==0
        u[99] += 1
      end
      if p[8]==0
        u[100] += 1
      end
  elseif temp<cum_rate[86]
      u[58] += 1
  elseif temp<cum_rate[87]
      u[57] += 1
  elseif temp<cum_rate[88]
      if maximum([p[7],p[8]])==0
        u[56] -= 1
      end
        u[102] -=1
  elseif temp<cum_rate[89]
      u[58] -= 1
  elseif temp<cum_rate[90]
      u[57] -= 1
  elseif temp<cum_rate[91]
      u[99] -= 1
  elseif temp<cum_rate[92]
      u[100] -= 1
  elseif temp<cum_rate[93]
      if maximum([p[7],p[8]])>1
        u[78] += 1
        u[56] -= 1
      else
        u[56] -=1
      end
      p[7]==1 ? u[99] += 1 : 0
      p[8]==1 ? u[100] += 1 : 0
  elseif temp<cum_rate[94]
      if maximum([p[7],p[8]])>2
        u[78] -= 1
        u[79] += 1
      else
        u[78] -= 1
      end
      p[7]==2 ? u[99] += 1 : 0
      p[8]==2 ? u[100] += 1 : 0
  elseif temp<cum_rate[95]
      if maximum([p[7],p[8]])>3
        u[79] -= 1
        u[80] += 1
      else
        u[79] -= 1
      end
      p[7]==3 ? u[99] += 1 : 0
      p[8]==3 ? u[100] += 1 : 0
  elseif temp<cum_rate[96]
      if maximum([p[7],p[8]])>4
        u[80] -= 1
        u[81] += 1
      else
        u[80] -= 1
      end
      p[7]==4 ? u[99] += 1 : 0
      p[8]==4 ? u[100] += 1 : 0
  elseif temp<cum_rate[97]
      if maximum([p[7],p[8]])>5
        u[81] -= 1
        u[82] += 1
      else
        u[81] -= 1
      end
      p[7]==5 ? u[99] += 1 : 0
      p[8]==5 ? u[100] += 1 : 0
  elseif temp<cum_rate[98]
      if maximum([p[7],p[8]])>6
        u[82] -= 1
        u[83] += 1
      else
        u[82] -= 1
      end
      p[7]==6 ? u[99] += 1 : 0
      p[8]==6 ? u[100] += 1 : 0
  elseif temp<cum_rate[99]
      if maximum([p[7],p[8]])>7
        u[83] -= 1
        u[84] += 1
      else
        u[83] -= 1
      end
      p[7]==7 ? u[99] += 1 : 0
      p[8]==7 ? u[100] += 1 : 0
  elseif temp<cum_rate[100]
      if maximum([p[7],p[8]])>8
        u[84] -= 1
        u[85] += 1
      else
        u[84] -= 1
      end
      p[7]==8 ? u[99] += 1 : 0
      p[8]==8 ? u[100] += 1 : 0
  elseif temp<cum_rate[101]
      if maximum([p[7],p[8]])>9
        u[85] -= 1
        u[86] += 1
      else
        u[85] -= 1
      end
      p[7]==9 ? u[99] += 1 : 0
      p[8]==9 ? u[100] += 1 : 0
  elseif temp<cum_rate[102]
      if maximum([p[7],p[8]])>10
        u[86] -= 1
        u[87] += 1
      else
        u[86] -= 1
      end
      p[7]==10 ? u[99] += 1 : 0
      p[8]==10 ? u[100] += 1 : 0
  elseif temp<cum_rate[103]
      if maximum([p[7],p[8]])>11
        u[87] -= 1
        u[88] += 1
      else
        u[87] -= 1
      end
      p[7]==11 ? u[99] += 1 : 0
      p[8]==11 ? u[100] += 1 : 0
  elseif temp<cum_rate[104]
      if maximum([p[7],p[8]])>12
        u[88] -= 1
        u[89] += 1
      else
        u[88] -= 1
      end
      p[7]==12 ? u[99] += 1 : 0
      p[8]==12 ? u[100] += 1 : 0
  elseif temp<cum_rate[105]
      if maximum([p[7],p[8]])>13
        u[89] -= 1
        u[90] += 1
      else
        u[89] -= 1
      end
      p[7]==13 ? u[99] += 1 : 0
      p[8]==13 ? u[100] += 1 : 0
  elseif temp<cum_rate[106]
      if maximum([p[7],p[8]])>14
        u[90] -= 1
        u[91] += 1
      else
        u[90] -= 1
      end
      p[7]==14 ? u[99] += 1 : 0
      p[8]==14 ? u[100] += 1 : 0
  elseif temp<cum_rate[107]
      if maximum([p[7],p[8]])>15
        u[91] -= 1
        u[92] += 1
      else
        u[91] -= 1
      end
      p[7]==15 ? u[99] += 1 : 0
      p[8]==15 ? u[100] += 1 : 0
  elseif temp<cum_rate[108]
      if maximum([p[7],p[8]])>16
        u[92] -= 1
        u[93] += 1
      else
        u[92] -= 1
      end
      p[7]==16 ? u[99] += 1 : 0
      p[8]==16 ? u[100] += 1 : 0
  elseif temp<cum_rate[109]
      if maximum([p[7],p[8]])>17
        u[93] -= 1
        u[94] += 1
      else
        u[93] -= 1
      end
      p[7]==17 ? u[99] += 1 : 0
      p[8]==17 ? u[100] += 1 : 0
  elseif temp<cum_rate[110]
      if maximum([p[7],p[8]])>18
        u[94] -= 1
        u[95] += 1
      else
        u[94] -= 1
      end
      p[7]==18 ? u[99] += 1 : 0
      p[8]==18 ? u[100] += 1 : 0
  elseif temp<cum_rate[111]
      if maximum([p[7],p[8]])>19
        u[95] -= 1
        u[96] += 1
      else
        u[95] -= 1
      end
      p[7]==19 ? u[99] += 1 : 0
      p[8]==19 ? u[100] += 1 : 0
  elseif temp<cum_rate[112]
      u[96] -= 1
      p[7]==20 ? u[99] += 1 : 0
      p[8]==20 ? u[100] += 1 : 0
  end
    t = t+t_step
    push!(out,copy(u))
    push!(t_out,copy(t))
  end
  return out, t_out
end


### Polyploid cells with XA dilution hypothesis
function sim_delay_polypl_dil(p::Array{Float64,1},tspan::Tuple{Float64,Float64},u0::Array{Int64,1})
  t::Float64=tspan[1]
  i::Int64=0
  out::Array{Array{Int64,1}} = []
  t_out::Array{Float64,1} = []
  push!(out, u0)
  push!(t_out,t)
  u::Array{Int64,1}=copy(u0)
  #############################################################
  rate::Array{Float64,1}=zeros(56)
  if p[26] == 3
      rate = [rate; zeros(28)]
  elseif p[26] ==4
      rate = [rate; zeros(56)]
  end
  p[24] = 0
  p[25] = 0
  p[27] = 0
  p[28] = 0
  p[29] = 1
  p[30] = 1
  p[31] = 1
  p[32] = 1
  while (t<tspan[2])
    i+=1
    if p[26] == 2
        # Xist transcription
        rate[1] = p[21]*(p[24]+p[29]*(1/p[33]*(u[5]+u[6]))^p[11]/((1/p[33]*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[7]^p[13]/(u[7]^p[13]+(p[22]*p[14])^p[13])))
        rate[2] = p[21]*(p[25]+p[30]*(1/p[33]*(u[5]+u[6]))^p[11]/((1/p[33]*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[8]^p[13]/(u[8]^p[13]+(p[22]*p[14])^p[13])))
        # cXR production
        rate[3] = p[22]*(1-u[51]^p[5]/(u[51]^p[5]+(p[21]*p[6])^p[5]))
        rate[4] = p[22]*(1-u[49]^p[5]/(u[49]^p[5]+(p[21]*p[6])^p[5]))
        # tXA production
        rate[5] = p[23]*(1-u[52]^p[3]/(u[52]^p[3]+(p[21]*p[4])^p[3]))
        rate[6] = p[23]*(1-u[50]^p[3]/(u[50]^p[3]+(p[21]*p[4])^p[3]))
        # degradation
        rate[7] = 0.1733*u[3]
        rate[8] = 0.1733*u[4]
        rate[9] = 1*u[7]
        rate[10] = 1*u[8]
        rate[11] = 1*u[5]
        rate[12] = 1*u[6]
        rate[13] = 1*u[49]
        rate[14] = 1*u[50]
        rate[15] = 1*u[51]
        rate[16] = 1*u[52]
        ### silencing intermediates on chr1
        rate[17] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[1],0)
        rate[18] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[10],0)
        rate[19] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[11],0)
        rate[20] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[12],0)
        rate[21] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[13],0)
        rate[22] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[14],0)
        rate[23] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[15],0)
        rate[24] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[16],0)
        rate[25] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[17],0)
        rate[26] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[18],0)
        rate[27] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[19],0)
        rate[28] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[20],0)
        rate[29] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[21],0)
        rate[30] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[22],0)
        rate[31] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[23],0)
        rate[32] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[24],0)
        rate[33] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[25],0)
        rate[34] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[26],0)
        rate[35] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[27],0)
        rate[36] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[28],0)
        ### silencing intermediates chr 2
        rate[37] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[2],0)
        rate[38] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[30],0)
        rate[39] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[31],0)
        rate[40] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[32],0)
        rate[41] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[33],0)
        rate[42] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[34],0)
        rate[43] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[35],0)
        rate[44] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[36],0)
        rate[45] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[37],0)
        rate[46] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[38],0)
        rate[47] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[39],0)
        rate[48] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[40],0)
        rate[49] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[41],0)
        rate[50] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[42],0)
        rate[51] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[43],0)
        rate[52] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[44],0)
        rate[53] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[45],0)
        rate[54] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[46],0)
        rate[55] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[47],0)
        rate[56] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[48],0)
    elseif p[26] == 3
        # Xist transcription
        rate[1] = p[21]*(p[24]+p[29]*(1/p[33]*(u[5]+u[6]+u[54]))^p[11]/((1/p[33]*(u[5]+u[6]+u[54]))^p[11]+(p[23]*p[12])^p[11])*(1-u[7]^p[13]/(u[7]^p[13]+(p[22]*p[14])^p[13])))
        rate[2] = p[21]*(p[25]+p[30]*(1/p[33]*(u[5]+u[6]+u[54]))^p[11]/((1/p[33]*(u[5]+u[6]+u[54]))^p[11]+(p[23]*p[12])^p[11])*(1-u[8]^p[13]/(u[8]^p[13]+(p[22]*p[14])^p[13])))
        rate[57] = p[21]*(p[27]+p[31]*(1/p[33]*(u[5]+u[6]+u[54]))^p[11]/((1/p[33]*(u[5]+u[6]+u[54]))^p[11]+(p[23]*p[12])^p[11])*(1-u[55]^p[13]/(u[55]^p[13]+(p[22]*p[14])^p[13])))
        # cXR production (scaling factor * (1-silencing))
        rate[3] = p[22]*(1-u[51]^p[5]/(u[51]^p[5]+(p[21]*p[6])^p[5]))
        rate[4] = p[22]*(1-u[49]^p[5]/(u[49]^p[5]+(p[21]*p[6])^p[5]))
        rate[58] = p[22]*(1-u[97]^p[5]/(u[97]^p[5]+(p[21]*p[6])^p[5]))
        # tXA production (scaling factor * (1-silencing))
        rate[5] = p[23]*(1-u[52]^p[3]/(u[52]^p[3]+(p[21]*p[4])^p[3]))
        rate[6] = p[23]*(1-u[50]^p[3]/(u[50]^p[3]+(p[21]*p[4])^p[3]))
        rate[59] = p[23]*(1-u[98]^p[3]/(u[98]^p[3]+(p[21]*p[4])^p[3]))
        # degradation
        rate[7] = 0.1733*u[3]
        rate[8] = 0.1733*u[4]
        rate[60] = 0.1733*u[101]
        rate[9] = 1*u[7]
        rate[10] = 1*u[8]
        rate[61] = 1*u[55]
        rate[11] = 1*u[5]
        rate[12] = 1*u[6]
        rate[62] = 1*u[54]
        rate[13] = 1*u[49]
        rate[14] = 1*u[50]
        rate[63] = 1*u[97]
        rate[15] = 1*u[51]
        rate[16] = 1*u[52]
        rate[64] = 1*u[98]
        ### silencing intermediates on chr1
        rate[17] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[1],0)
        rate[18] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[10],0)
        rate[19] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[11],0)
        rate[20] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[12],0)
        rate[21] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[13],0)
        rate[22] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[14],0)
        rate[23] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[15],0)
        rate[24] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[16],0)
        rate[25] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[17],0)
        rate[26] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[18],0)
        rate[27] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[19],0)
        rate[28] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[20],0)
        rate[29] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[21],0)
        rate[30] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[22],0)
        rate[31] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[23],0)
        rate[32] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[24],0)
        rate[33] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[25],0)
        rate[34] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[26],0)
        rate[35] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[27],0)
        rate[36] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[28],0)
        ### silencing intermediates chr 2
        rate[37] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[2],0)
        rate[38] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[30],0)
        rate[39] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[31],0)
        rate[40] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[32],0)
        rate[41] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[33],0)
        rate[42] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[34],0)
        rate[43] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[35],0)
        rate[44] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[36],0)
        rate[45] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[37],0)
        rate[46] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[38],0)
        rate[47] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[39],0)
        rate[48] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[40],0)
        rate[49] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[41],0)
        rate[50] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[42],0)
        rate[51] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[43],0)
        rate[52] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[44],0)
        rate[53] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[45],0)
        rate[54] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[46],0)
        rate[55] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[47],0)
        rate[56] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[48],0)
        ### silencing intermediates chr 3
        rate[65] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[53],0)
        rate[66] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[59],0)
        rate[67] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[60],0)
        rate[68] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[61],0)
        rate[69] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[62],0)
        rate[70] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[63],0)
        rate[71] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[64],0)
        rate[72] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[65],0)
        rate[73] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[66],0)
        rate[74] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[67],0)
        rate[75] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[68],0)
        rate[76] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[69],0)
        rate[77] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[70],0)
        rate[78] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[71],0)
        rate[79] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[72],0)
        rate[80] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[73],0)
        rate[81] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[74],0)
        rate[82] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[75],0)
        rate[83] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[76],0)
        rate[84] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[77],0)
    elseif p[26] == 4
        # Xist transcription
        rate[1] = p[21]*(p[24]+p[29]*(1/p[33]*(u[5]+u[6]+u[54]+u[57]))^p[11]/((1/p[33]*(u[5]+u[6]+u[54]+u[57]))^p[11]+(p[23]*p[12])^p[11])*(1-u[7]^p[13]/(u[7]^p[13]+(p[22]*p[14])^p[13])))
        rate[2] = p[21]*(p[25]+p[30]*(1/p[33]*(u[5]+u[6]+u[54]+u[57]))^p[11]/((1/p[33]*(u[5]+u[6]+u[54]+u[57]))^p[11]+(p[23]*p[12])^p[11])*(1-u[8]^p[13]/(u[8]^p[13]+(p[22]*p[14])^p[13])))
        rate[57] = p[21]*(p[27]+p[31]*(1/p[33]*(u[5]+u[6]+u[54]+u[57]))^p[11]/((1/p[33]*(u[5]+u[6]+u[54]+u[57]))^p[11]+(p[23]*p[12])^p[11])*(1-u[55]^p[13]/(u[55]^p[13]+(p[22]*p[14])^p[13])))
        rate[85] = p[21]*(p[28]+p[32]*(1/p[33]*(u[5]+u[6]+u[54]+u[57]))^p[11]/((1/p[33]*(u[5]+u[6]+u[54]+u[57]))^p[11]+(p[23]*p[12])^p[11])*(1-u[58]^p[13]/(u[58]^p[13]+(p[22]*p[14])^p[13])))
        # cXR production (scaling factor * (1-silencing))
        rate[3] = p[22]*(1-u[51]^p[5]/(u[51]^p[5]+(p[21]*p[6])^p[5]))
        rate[4] = p[22]*(1-u[49]^p[5]/(u[49]^p[5]+(p[21]*p[6])^p[5]))
        rate[58] = p[22]*(1-u[97]^p[5]/(u[97]^p[5]+(p[21]*p[6])^p[5]))
        rate[86] = p[22]*(1-u[99]^p[5]/(u[99]^p[5]+(p[21]*p[6])^p[5]))
        # tXA production (scaling factor * (1-silencing))
        rate[5] = p[23]*(1-u[52]^p[3]/(u[52]^p[3]+(p[21]*p[4])^p[3]))
        rate[6] = p[23]*(1-u[50]^p[3]/(u[50]^p[3]+(p[21]*p[4])^p[3]))
        rate[59] = p[23]*(1-u[98]^p[3]/(u[98]^p[3]+(p[21]*p[4])^p[3]))
        rate[87] = p[23]*(1-u[100]^p[3]/(u[100]^p[3]+(p[21]*p[4])^p[3]))
        # degradation
        rate[7] = 0.1733*u[3]
        rate[8] = 0.1733*u[4]
        rate[60] = 0.1733*u[101]
        rate[88] = 0.1733*u[102]
        rate[9] = 1*u[7]
        rate[10] = 1*u[8]
        rate[61] = 1*u[55]
        rate[89] = 1*u[58]
        rate[11] = 1*u[5]
        rate[12] = 1*u[6]
        rate[62] = 1*u[54]
        rate[90] = 1*u[57]
        rate[13] = 1*u[49]
        rate[14] = 1*u[50]
        rate[63] = 1*u[97]
        rate[91] = 1*u[99]
        rate[15] = 1*u[51]
        rate[16] = 1*u[52]
        rate[64] = 1*u[98]
        rate[92] = 1*u[100]
        ### silencing intermediates on chr1
        rate[17] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[1],0)
        rate[18] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[10],0)
        rate[19] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[11],0)
        rate[20] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[12],0)
        rate[21] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[13],0)
        rate[22] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[14],0)
        rate[23] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[15],0)
        rate[24] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[16],0)
        rate[25] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[17],0)
        rate[26] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[18],0)
        rate[27] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[19],0)
        rate[28] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[20],0)
        rate[29] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[21],0)
        rate[30] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[22],0)
        rate[31] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[23],0)
        rate[32] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[24],0)
        rate[33] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[25],0)
        rate[34] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[26],0)
        rate[35] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[27],0)
        rate[36] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[28],0)
        ### silencing intermediates chr 2
        rate[37] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[2],0)
        rate[38] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[30],0)
        rate[39] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[31],0)
        rate[40] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[32],0)
        rate[41] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[33],0)
        rate[42] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[34],0)
        rate[43] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[35],0)
        rate[44] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[36],0)
        rate[45] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[37],0)
        rate[46] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[38],0)
        rate[47] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[39],0)
        rate[48] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[40],0)
        rate[49] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[41],0)
        rate[50] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[42],0)
        rate[51] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[43],0)
        rate[52] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[44],0)
        rate[53] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[45],0)
        rate[54] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[46],0)
        rate[55] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[47],0)
        rate[56] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[48],0)
        ### silencing intermediates chr 3
        rate[65] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[53],0)
        rate[66] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[59],0)
        rate[67] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[60],0)
        rate[68] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[61],0)
        rate[69] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[62],0)
        rate[70] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[63],0)
        rate[71] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[64],0)
        rate[72] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[65],0)
        rate[73] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[66],0)
        rate[74] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[67],0)
        rate[75] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[68],0)
        rate[76] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[69],0)
        rate[77] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[70],0)
        rate[78] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[71],0)
        rate[79] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[72],0)
        rate[80] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[73],0)
        rate[81] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[74],0)
        rate[82] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[75],0)
        rate[83] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[76],0)
        rate[84] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[77],0)
        ### silencing intermediates chr 4
        rate[93] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[56],0)
        rate[94] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[78],0)
        rate[95] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[79],0)
        rate[96] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[80],0)
        rate[97] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[81],0)
        rate[98] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[82],0)
        rate[99] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[83],0)
        rate[100] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[84],0)
        rate[101] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[85],0)
        rate[102] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[86],0)
        rate[103] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[87],0)
        rate[104] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[88],0)
        rate[105] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[89],0)
        rate[106] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[90],0)
        rate[107] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[91],0)
        rate[108] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[92],0)
        rate[109] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[93],0)
        rate[110] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[94],0)
        rate[111] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[95],0)
        rate[112] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[96],0)
    end

    # Gillespie
    r = rand(2)
    # determine time point of next reaction
    t_step::Float64 = 1/sum(rate)*log(1/(1-r[1]))
    # decide which reaction will occur
    temp::Float64 = r[2]*sum(rate)
    cum_rate = cumsum(rate)
    # Execute reaction
    if temp<cum_rate[1]
      u[1] += 1
      u[3] += 1
        if p[7]==0
        u[49] += 1
      end
      if p[8]==0
        u[50] += 1
      end
  elseif temp<cum_rate[2]
      u[2] += 1
      u[4] += 1
      if p[7]==0
        u[51] += 1
      end
      if p[8]==0
        u[52] += 1
      end
  elseif temp<cum_rate[3]
        u[8] += 1
  elseif temp<cum_rate[4]
      u[7] += 1
  elseif temp<cum_rate[5]
        u[6] += 1
  elseif temp<cum_rate[6]
        u[5] += 1
  elseif temp<cum_rate[7]
      if maximum([p[7],p[8]])==0
        u[1] -= 1
      end
        u[3] -=1
  elseif temp<cum_rate[8]
      if maximum([p[7],p[8]])==0
        u[2] -= 1
      end
        u[4] -=1
  elseif temp<cum_rate[9]
        u[7] -= 1
  elseif temp<cum_rate[10]
        u[8] -= 1
  elseif temp<cum_rate[11]
        u[5] -= 1
  elseif temp<cum_rate[12]
        u[6] -= 1
  elseif temp<cum_rate[13]
        u[49] -= 1
  elseif temp<cum_rate[14]
        u[50] -= 1
  elseif temp<cum_rate[15]
        u[51] -= 1
  elseif temp<cum_rate[16]
        u[52] -= 1
  elseif temp<cum_rate[17]
      if maximum([p[7],p[8]])>1
        u[10] += 1
        u[1] -= 1
      else
        u[1] -= 1
      end
      p[7]==1 ? u[49] += 1 : 0
      p[8]==1 ? u[50] += 1 : 0
  elseif temp<cum_rate[18]
      if maximum([p[7],p[8]])>2
        u[10] -= 1
        u[11] += 1
      else
        u[10] -= 1
      end
      p[7]==2 ? u[49] += 1 : 0
      p[8]==2 ? u[50] += 1 : 0
  elseif temp<cum_rate[19]
      if maximum([p[7],p[8]])>3
        u[11] -= 1
        u[12] += 1
      else
        u[11] -= 1
      end
      p[7]==3 ? u[49] += 1 : 0
      p[8]==3 ? u[50] += 1 : 0
  elseif temp<cum_rate[20]
      if maximum([p[7],p[8]])>4
        u[12] -= 1
        u[13] += 1
      else
        u[12] -= 1
      end
      p[7]==4 ? u[49] += 1 : 0
      p[8]==4 ? u[50] += 1 : 0
  elseif temp<cum_rate[21]
      if maximum([p[7],p[8]])>5
        u[13] -= 1
        u[14] += 1
      else
        u[13] -= 1
      end
      p[7]==5 ? u[49] += 1 : 0
      p[8]==5 ? u[50] += 1 : 0
  elseif temp<cum_rate[22]
      if maximum([p[7],p[8]])>6
        u[14] -= 1
        u[15] += 1
      else
        u[14] -= 1
      end
      p[7]==6 ? u[49] += 1 : 0
      p[8]==6 ? u[50] += 1 : 0
  elseif temp<cum_rate[23]
      if maximum([p[7],p[8]])>7
        u[15] -= 1
        u[16] += 1
      else
        u[15] -= 1
      end
      p[7]==7 ? u[49] += 1 : 0
      p[8]==7 ? u[50] += 1 : 0
  elseif temp<cum_rate[24]
      if maximum([p[7],p[8]])>8
        u[16] -= 1
        u[17] += 1
      else
        u[16] -= 1
      end
      p[7]==8 ? u[49] += 1 : 0
      p[8]==8 ? u[50] += 1 : 0
  elseif temp<cum_rate[25]
      if maximum([p[7],p[8]])>9
        u[17] -= 1
        u[18] += 1
      else
        u[17] -= 1
      end
      p[7]==9 ? u[49] += 1 : 0
      p[8]==9 ? u[50] += 1 : 0
  elseif temp<cum_rate[26]
      if maximum([p[7],p[8]])>10
        u[18] -= 1
        u[19] += 1
      else
        u[18] -= 1
      end
      p[7]==10 ? u[49] += 1 : 0
      p[8]==10 ? u[50] += 1 : 0
  elseif temp<cum_rate[27]
      if maximum([p[7],p[8]])>11
        u[19] -= 1
        u[20] += 1
      else
        u[19] -= 1
      end
      p[7]==11 ? u[49] += 1 : 0
      p[8]==11 ? u[50] += 1 : 0
  elseif temp<cum_rate[28]
      if maximum([p[7],p[8]])>12
        u[20] -= 1
        u[21] += 1
      else
        u[20] -= 1
      end
      p[7]==12 ? u[49] += 1 : 0
      p[8]==12 ? u[50] += 1 : 0
  elseif temp<cum_rate[29]
      if maximum([p[7],p[8]])>13
        u[21] -= 1
        u[22] += 1
      else
        u[21] -= 1
      end
      p[7]==13 ? u[49] += 1 : 0
      p[8]==13 ? u[50] += 1 : 0
  elseif temp<cum_rate[30]
      if maximum([p[7],p[8]])>14
        u[22] -= 1
        u[23] += 1
      else
        u[22] -= 1
      end
      p[7]==14 ? u[49] += 1 : 0
      p[8]==14 ? u[50] += 1 : 0
  elseif temp<cum_rate[31]
      if maximum([p[7],p[8]])>15
        u[23] -= 1
        u[24] += 1
      else
        u[23] -= 1
      end
      p[7]==15 ? u[49] += 1 : 0
      p[8]==15 ? u[50] += 1 : 0
  elseif temp<cum_rate[32]
      if maximum([p[7],p[8]])>16
        u[24] -= 1
        u[25] += 1
      else
        u[24] -= 1
      end
      p[7]==16 ? u[49] += 1 : 0
      p[8]==16 ? u[50] += 1 : 0
  elseif temp<cum_rate[33]
      if maximum([p[7],p[8]])>17
        u[25] -= 1
        u[26] += 1
      else
        u[25] -= 1
      end
      p[7]==17 ? u[49] += 1 : 0
      p[8]==17 ? u[50] += 1 : 0
  elseif temp<cum_rate[34]
      if maximum([p[7],p[8]])>18
        u[26] -= 1
        u[27] += 1
      else
        u[26] -= 1
      end
      p[7]==18 ? u[49] += 1 : 0
      p[8]==18 ? u[50] += 1 : 0
  elseif temp<cum_rate[35]
      if maximum([p[7],p[8]])>19
        u[27] -= 1
        u[28] += 1
      else
        u[27] -= 1
      end
      p[7]==19 ? u[49] += 1 : 0
      p[8]==19 ? u[50] += 1 : 0
  elseif temp<cum_rate[36]
      u[28] -= 1
      p[7]==20 ? u[49] += 1 : 0
      p[8]==20 ? u[50] += 1 : 0
  elseif temp<cum_rate[37]
      if maximum([p[7],p[8]])>1
        u[30] += 1
        u[2] -= 1
      else
        u[2] -=1
      end
      p[7]==1 ? u[51] += 1 : 0
      p[8]==1 ? u[52] += 1 : 0
  elseif temp<cum_rate[38]
      if maximum([p[7],p[8]])>2
        u[30] -= 1
        u[31] += 1
      else
        u[30] -= 1
      end
      p[7]==2 ? u[51] += 1 : 0
      p[8]==2 ? u[52] += 1 : 0
  elseif temp<cum_rate[39]
      if maximum([p[7],p[8]])>3
        u[31] -= 1
        u[32] += 1
      else
        u[31] -= 1
      end
      p[7]==3 ? u[51] += 1 : 0
      p[8]==3 ? u[52] += 1 : 0
  elseif temp<cum_rate[40]
      if maximum([p[7],p[8]])>4
        u[32] -= 1
        u[33] += 1
      else
        u[32] -= 1
      end
      p[7]==4 ? u[51] += 1 : 0
      p[8]==4 ? u[52] += 1 : 0
  elseif temp<cum_rate[41]
      if maximum([p[7],p[8]])>5
        u[33] -= 1
        u[34] += 1
      else
        u[33] -= 1
      end
      p[7]==5 ? u[51] += 1 : 0
      p[8]==5 ? u[52] += 1 : 0
  elseif temp<cum_rate[42]
      if maximum([p[7],p[8]])>6
        u[34] -= 1
        u[35] += 1
      else
        u[34] -= 1
      end
      p[7]==6 ? u[51] += 1 : 0
      p[8]==6 ? u[52] += 1 : 0
  elseif temp<cum_rate[43]
      if maximum([p[7],p[8]])>7
        u[35] -= 1
        u[36] += 1
      else
        u[35] -= 1
      end
      p[7]==7 ? u[51] += 1 : 0
      p[8]==7 ? u[52] += 1 : 0
  elseif temp<cum_rate[44]
      if maximum([p[7],p[8]])>8
        u[36] -= 1
        u[37] += 1
      else
        u[36] -= 1
      end
      p[7]==8 ? u[51] += 1 : 0
      p[8]==8 ? u[52] += 1 : 0
  elseif temp<cum_rate[45]
      if maximum([p[7],p[8]])>9
        u[37] -= 1
        u[38] += 1
      else
        u[37] -= 1
      end
      p[7]==9 ? u[51] += 1 : 0
      p[8]==9 ? u[52] += 1 : 0
  elseif temp<cum_rate[46]
      if maximum([p[7],p[8]])>10
        u[38] -= 1
        u[39] += 1
      else
        u[38] -= 1
      end
      p[7]==10 ? u[51] += 1 : 0
      p[8]==10 ? u[52] += 1 : 0
  elseif temp<cum_rate[47]
      if maximum([p[7],p[8]])>11
        u[39] -= 1
        u[40] += 1
      else
        u[39] -= 1
      end
      p[7]==11 ? u[51] += 1 : 0
      p[8]==11 ? u[52] += 1 : 0
  elseif temp<cum_rate[48]
      if maximum([p[7],p[8]])>12
        u[40] -= 1
        u[41] += 1
      else
        u[40] -= 1
      end
      p[7]==12 ? u[51] += 1 : 0
      p[8]==12 ? u[52] += 1 : 0
  elseif temp<cum_rate[49]
      if maximum([p[7],p[8]])>13
        u[41] -= 1
        u[42] += 1
      else
        u[41] -= 1
      end
      p[7]==13 ? u[51] += 1 : 0
      p[8]==13 ? u[52] += 1 : 0
  elseif temp<cum_rate[50]
      if maximum([p[7],p[8]])>14
        u[42] -= 1
        u[43] += 1
      else
        u[42] -= 1
      end
      p[7]==14 ? u[51] += 1 : 0
      p[8]==14 ? u[52] += 1 : 0
  elseif temp<cum_rate[51]
      if maximum([p[7],p[8]])>15
        u[43] -= 1
        u[44] += 1
      else
        u[43] -= 1
      end
      p[7]==15 ? u[51] += 1 : 0
      p[8]==15 ? u[52] += 1 : 0
  elseif temp<cum_rate[52]
      if maximum([p[7],p[8]])>16
        u[44] -= 1
        u[45] += 1
      else
        u[44] -= 1
      end
      p[7]==16 ? u[51] += 1 : 0
      p[8]==16 ? u[52] += 1 : 0
  elseif temp<cum_rate[53]
      if maximum([p[7],p[8]])>17
        u[45] -= 1
        u[46] += 1
      else
        u[45] -= 1
      end
      p[7]==17 ? u[51] += 1 : 0
      p[8]==17 ? u[52] += 1 : 0
  elseif temp<cum_rate[54]
      if maximum([p[7],p[8]])>18
        u[46] -= 1
        u[47] += 1
      else
        u[46] -= 1
      end
      p[7]==18 ? u[51] += 1 : 0
      p[8]==18 ? u[52] += 1 : 0
  elseif temp<cum_rate[55]
      if maximum([p[7],p[8]])>19
        u[47] -= 1
        u[48] += 1
      else
        u[47] -= 1
      end
      p[7]==19 ? u[51] += 1 : 0
      p[8]==19 ? u[52] += 1 : 0
  elseif temp<cum_rate[56]
      u[48] -= 1
      p[7]==20 ? u[51] += 1 : 0
      p[8]==20 ? u[52] += 1 : 0
  elseif temp<cum_rate[57]
      u[53] += 1
      u[101] += 1
      if p[7]==0
        u[97] += 1
      end
      if p[8]==0
        u[98] += 1
      end
  elseif temp<cum_rate[58]
      u[55] += 1
  elseif temp<cum_rate[59]
      u[54] += 1
  elseif temp<cum_rate[60]
      if maximum([p[7],p[8]])==0
        u[53] -= 1
      end
        u[101] -=1
  elseif temp<cum_rate[61]
      u[55] -= 1
  elseif temp<cum_rate[62]
      u[54] -= 1
  elseif temp<cum_rate[63]
      u[97] -= 1
  elseif temp<cum_rate[64]
      u[98] -= 1
  elseif temp<cum_rate[65]
      if maximum([p[7],p[8]])>1
        u[59] += 1
        u[53] -= 1
      else
        u[53] -=1
      end
      p[7]==1 ? u[97] += 1 : 0
      p[8]==1 ? u[98] += 1 : 0
  elseif temp<cum_rate[66]
      if maximum([p[7],p[8]])>2
        u[59] -= 1
        u[60] += 1
      else
        u[59] -= 1
      end
      p[7]==2 ? u[97] += 1 : 0
      p[8]==2 ? u[98] += 1 : 0
  elseif temp<cum_rate[67]
      if maximum([p[7],p[8]])>3
        u[60] -= 1
        u[61] += 1
      else
        u[60] -= 1
      end
      p[7]==3 ? u[97] += 1 : 0
      p[8]==3 ? u[98] += 1 : 0
  elseif temp<cum_rate[68]
      if maximum([p[7],p[8]])>4
        u[61] -= 1
        u[62] += 1
      else
        u[61] -= 1
      end
      p[7]==4 ? u[97] += 1 : 0
      p[8]==4 ? u[98] += 1 : 0
  elseif temp<cum_rate[69]
      if maximum([p[7],p[8]])>5
        u[62] -= 1
        u[63] += 1
      else
        u[62] -= 1
      end
      p[7]==5 ? u[97] += 1 : 0
      p[8]==5 ? u[98] += 1 : 0
  elseif temp<cum_rate[70]
      if maximum([p[7],p[8]])>6
        u[63] -= 1
        u[64] += 1
      else
        u[63] -= 1
      end
      p[7]==6 ? u[97] += 1 : 0
      p[8]==6 ? u[98] += 1 : 0
  elseif temp<cum_rate[71]
      if maximum([p[7],p[8]])>7
        u[64] -= 1
        u[65] += 1
      else
        u[64] -= 1
      end
      p[7]==7 ? u[97] += 1 : 0
      p[8]==7 ? u[98] += 1 : 0
  elseif temp<cum_rate[72]
      if maximum([p[7],p[8]])>8
        u[65] -= 1
        u[66] += 1
      else
        u[65] -= 1
      end
      p[7]==8 ? u[97] += 1 : 0
      p[8]==8 ? u[98] += 1 : 0
  elseif temp<cum_rate[73]
      if maximum([p[7],p[8]])>9
        u[66] -= 1
        u[67] += 1
      else
        u[66] -= 1
      end
      p[7]==9 ? u[97] += 1 : 0
      p[8]==9 ? u[98] += 1 : 0
  elseif temp<cum_rate[74]
      if maximum([p[7],p[8]])>10
        u[67] -= 1
        u[68] += 1
      else
        u[67] -= 1
      end
      p[7]==10 ? u[97] += 1 : 0
      p[8]==10 ? u[98] += 1 : 0
  elseif temp<cum_rate[75]
      if maximum([p[7],p[8]])>11
        u[68] -= 1
        u[69] += 1
      else
        u[68] -= 1
      end
      p[7]==11 ? u[97] += 1 : 0
      p[8]==11 ? u[98] += 1 : 0
  elseif temp<cum_rate[76]
      if maximum([p[7],p[8]])>12
        u[69] -= 1
        u[70] += 1
      else
        u[69] -= 1
      end
      p[7]==12 ? u[97] += 1 : 0
      p[8]==12 ? u[98] += 1 : 0
  elseif temp<cum_rate[77]
      if maximum([p[7],p[8]])>13
        u[70] -= 1
        u[71] += 1
      else
        u[70] -= 1
      end
      p[7]==13 ? u[97] += 1 : 0
      p[8]==13 ? u[98] += 1 : 0
  elseif temp<cum_rate[78]
      if maximum([p[7],p[8]])>14
        u[71] -= 1
        u[72] += 1
      else
        u[71] -= 1
      end
      p[7]==14 ? u[97] += 1 : 0
      p[8]==14 ? u[98] += 1 : 0
  elseif temp<cum_rate[79]
      if maximum([p[7],p[8]])>15
        u[72] -= 1
        u[73] += 1
      else
        u[72] -= 1
      end
      p[7]==15 ? u[97] += 1 : 0
      p[8]==15 ? u[98] += 1 : 0
  elseif temp<cum_rate[80]
      if maximum([p[7],p[8]])>16
        u[73] -= 1
        u[74] += 1
      else
        u[73] -= 1
      end
      p[7]==16 ? u[97] += 1 : 0
      p[8]==16 ? u[98] += 1 : 0
  elseif temp<cum_rate[81]
      if maximum([p[7],p[8]])>17
        u[74] -= 1
        u[75] += 1
      else
        u[74] -= 1
      end
      p[7]==17 ? u[97] += 1 : 0
      p[8]==17 ? u[98] += 1 : 0
  elseif temp<cum_rate[82]
      if maximum([p[7],p[8]])>18
        u[75] -= 1
        u[76] += 1
      else
        u[75] -= 1
      end
      p[7]==18 ? u[97] += 1 : 0
      p[8]==18 ? u[98] += 1 : 0
  elseif temp<cum_rate[83]
      if maximum([p[7],p[8]])>19
        u[76] -= 1
        u[77] += 1
      else
        u[76] -= 1
      end
      p[7]==19 ? u[97] += 1 : 0
      p[8]==19 ? u[98] += 1 : 0
  elseif temp<cum_rate[84]
      u[77] -= 1
      p[7]==20 ? u[97] += 1 : 0
      p[8]==20 ? u[98] += 1 : 0
  elseif temp<cum_rate[85]
      u[56] += 1
      u[102] += 1
      if p[7]==0
        u[99] += 1
      end
      if p[8]==0
        u[100] += 1
      end
  elseif temp<cum_rate[86]
      u[58] += 1
  elseif temp<cum_rate[87]
      u[57] += 1
  elseif temp<cum_rate[88]
      if maximum([p[7],p[8]])==0
        u[56] -= 1
      end
        u[102] -=1
  elseif temp<cum_rate[89]
      u[58] -= 1
  elseif temp<cum_rate[90]
      u[57] -= 1
  elseif temp<cum_rate[91]
      u[99] -= 1
  elseif temp<cum_rate[92]
      u[100] -= 1
  elseif temp<cum_rate[93]
      if maximum([p[7],p[8]])>1
        u[78] += 1
        u[56] -= 1
      else
        u[56] -=1
      end
      p[7]==1 ? u[99] += 1 : 0
      p[8]==1 ? u[100] += 1 : 0
  elseif temp<cum_rate[94]
      if maximum([p[7],p[8]])>2
        u[78] -= 1
        u[79] += 1
      else
        u[78] -= 1
      end
      p[7]==2 ? u[99] += 1 : 0
      p[8]==2 ? u[100] += 1 : 0
  elseif temp<cum_rate[95]
      if maximum([p[7],p[8]])>3
        u[79] -= 1
        u[80] += 1
      else
        u[79] -= 1
      end
      p[7]==3 ? u[99] += 1 : 0
      p[8]==3 ? u[100] += 1 : 0
  elseif temp<cum_rate[96]
      if maximum([p[7],p[8]])>4
        u[80] -= 1
        u[81] += 1
      else
        u[80] -= 1
      end
      p[7]==4 ? u[99] += 1 : 0
      p[8]==4 ? u[100] += 1 : 0
  elseif temp<cum_rate[97]
      if maximum([p[7],p[8]])>5
        u[81] -= 1
        u[82] += 1
      else
        u[81] -= 1
      end
      p[7]==5 ? u[99] += 1 : 0
      p[8]==5 ? u[100] += 1 : 0
  elseif temp<cum_rate[98]
      if maximum([p[7],p[8]])>6
        u[82] -= 1
        u[83] += 1
      else
        u[82] -= 1
      end
      p[7]==6 ? u[99] += 1 : 0
      p[8]==6 ? u[100] += 1 : 0
  elseif temp<cum_rate[99]
      if maximum([p[7],p[8]])>7
        u[83] -= 1
        u[84] += 1
      else
        u[83] -= 1
      end
      p[7]==7 ? u[99] += 1 : 0
      p[8]==7 ? u[100] += 1 : 0
  elseif temp<cum_rate[100]
      if maximum([p[7],p[8]])>8
        u[84] -= 1
        u[85] += 1
      else
        u[84] -= 1
      end
      p[7]==8 ? u[99] += 1 : 0
      p[8]==8 ? u[100] += 1 : 0
  elseif temp<cum_rate[101]
      if maximum([p[7],p[8]])>9
        u[85] -= 1
        u[86] += 1
      else
        u[85] -= 1
      end
      p[7]==9 ? u[99] += 1 : 0
      p[8]==9 ? u[100] += 1 : 0
  elseif temp<cum_rate[102]
      if maximum([p[7],p[8]])>10
        u[86] -= 1
        u[87] += 1
      else
        u[86] -= 1
      end
      p[7]==10 ? u[99] += 1 : 0
      p[8]==10 ? u[100] += 1 : 0
  elseif temp<cum_rate[103]
      if maximum([p[7],p[8]])>11
        u[87] -= 1
        u[88] += 1
      else
        u[87] -= 1
      end
      p[7]==11 ? u[99] += 1 : 0
      p[8]==11 ? u[100] += 1 : 0
  elseif temp<cum_rate[104]
      if maximum([p[7],p[8]])>12
        u[88] -= 1
        u[89] += 1
      else
        u[88] -= 1
      end
      p[7]==12 ? u[99] += 1 : 0
      p[8]==12 ? u[100] += 1 : 0
  elseif temp<cum_rate[105]
      if maximum([p[7],p[8]])>13
        u[89] -= 1
        u[90] += 1
      else
        u[89] -= 1
      end
      p[7]==13 ? u[99] += 1 : 0
      p[8]==13 ? u[100] += 1 : 0
  elseif temp<cum_rate[106]
      if maximum([p[7],p[8]])>14
        u[90] -= 1
        u[91] += 1
      else
        u[90] -= 1
      end
      p[7]==14 ? u[99] += 1 : 0
      p[8]==14 ? u[100] += 1 : 0
  elseif temp<cum_rate[107]
      if maximum([p[7],p[8]])>15
        u[91] -= 1
        u[92] += 1
      else
        u[91] -= 1
      end
      p[7]==15 ? u[99] += 1 : 0
      p[8]==15 ? u[100] += 1 : 0
  elseif temp<cum_rate[108]
      if maximum([p[7],p[8]])>16
        u[92] -= 1
        u[93] += 1
      else
        u[92] -= 1
      end
      p[7]==16 ? u[99] += 1 : 0
      p[8]==16 ? u[100] += 1 : 0
  elseif temp<cum_rate[109]
      if maximum([p[7],p[8]])>17
        u[93] -= 1
        u[94] += 1
      else
        u[93] -= 1
      end
      p[7]==17 ? u[99] += 1 : 0
      p[8]==17 ? u[100] += 1 : 0
  elseif temp<cum_rate[110]
      if maximum([p[7],p[8]])>18
        u[94] -= 1
        u[95] += 1
      else
        u[94] -= 1
      end
      p[7]==18 ? u[99] += 1 : 0
      p[8]==18 ? u[100] += 1 : 0
  elseif temp<cum_rate[111]
      if maximum([p[7],p[8]])>19
        u[95] -= 1
        u[96] += 1
      else
        u[95] -= 1
      end
      p[7]==19 ? u[99] += 1 : 0
      p[8]==19 ? u[100] += 1 : 0
  elseif temp<cum_rate[112]
      u[96] -= 1
      p[7]==20 ? u[99] += 1 : 0
      p[8]==20 ? u[100] += 1 : 0
  end
    t = t+t_step
    push!(out,copy(u))
    push!(t_out,copy(t))
  end
  return out, t_out
end

### Polyploid cells with XA repression hypothesis
function sim_delay_polypl_repXA(p::Array{Float64,1},tspan::Tuple{Float64,Float64},u0::Array{Int64,1})
  t::Float64=tspan[1]
  i::Int64=0
  out::Array{Array{Int64,1}} = []
  t_out::Array{Float64,1} = []
  push!(out, u0)
  push!(t_out,t)
  u::Array{Int64,1}=copy(u0)
  #############################################################
  rate::Array{Float64,1}=zeros(56)
  if p[26] == 3
      rate = [rate; zeros(28)]
  elseif p[26] ==4
      rate = [rate; zeros(56)]
  end
  p[24] = 0
  p[25] = 0
  p[27] = 0
  p[28] = 0
  p[29] = 1
  p[30] = 1
  p[31] = 1
  p[32] = 1
  while (t<tspan[2])
    i+=1
    if p[26] == 2
        # Xist transcription
        rate[1] = p[21]*(p[24]+p[29]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[7]^p[13]/(u[7]^p[13]+(p[22]*p[14])^p[13])))
        rate[2] = p[21]*(p[25]+p[30]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[8]^p[13]/(u[8]^p[13]+(p[22]*p[14])^p[13])))
        # cXR production
        rate[3] = p[22]*(1-u[51]^p[5]/(u[51]^p[5]+(p[21]*p[6])^p[5]))
        rate[4] = p[22]*(1-u[49]^p[5]/(u[49]^p[5]+(p[21]*p[6])^p[5]))
        # tXA production
        rate[5] = p[23]*2/p[33]*(1-u[52]^p[3]/(u[52]^p[3]+(p[21]*p[4])^p[3]))
        rate[6] = p[23]*2/p[33]*(1-u[50]^p[3]/(u[50]^p[3]+(p[21]*p[4])^p[3]))
        # degradation
        rate[7] = 0.1733*u[3]
        rate[8] = 0.1733*u[4]
        rate[9] = 1*u[7]
        rate[10] = 1*u[8]
        rate[11] = 1*u[5]
        rate[12] = 1*u[6]
        rate[13] = 1*u[49]
        rate[14] = 1*u[50]
        rate[15] = 1*u[51]
        rate[16] = 1*u[52]
        ### silencing intermediates on chr1
        rate[17] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[1],0)
        rate[18] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[10],0)
        rate[19] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[11],0)
        rate[20] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[12],0)
        rate[21] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[13],0)
        rate[22] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[14],0)
        rate[23] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[15],0)
        rate[24] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[16],0)
        rate[25] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[17],0)
        rate[26] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[18],0)
        rate[27] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[19],0)
        rate[28] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[20],0)
        rate[29] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[21],0)
        rate[30] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[22],0)
        rate[31] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[23],0)
        rate[32] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[24],0)
        rate[33] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[25],0)
        rate[34] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[26],0)
        rate[35] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[27],0)
        rate[36] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[28],0)
        ### silencing intermediates chr 2
        rate[37] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[2],0)
        rate[38] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[30],0)
        rate[39] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[31],0)
        rate[40] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[32],0)
        rate[41] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[33],0)
        rate[42] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[34],0)
        rate[43] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[35],0)
        rate[44] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[36],0)
        rate[45] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[37],0)
        rate[46] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[38],0)
        rate[47] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[39],0)
        rate[48] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[40],0)
        rate[49] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[41],0)
        rate[50] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[42],0)
        rate[51] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[43],0)
        rate[52] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[44],0)
        rate[53] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[45],0)
        rate[54] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[46],0)
        rate[55] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[47],0)
        rate[56] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[48],0)
    elseif p[26] == 3
        # Xist transcription
        rate[1] = p[21]*(p[24]+p[29]*(0.5*(u[5]+u[6]+u[54]))^p[11]/((0.5*(u[5]+u[6]+u[54]))^p[11]+(p[23]*p[12])^p[11])*(1-u[7]^p[13]/(u[7]^p[13]+(p[22]*p[14])^p[13])))
        rate[2] = p[21]*(p[25]+p[30]*(0.5*(u[5]+u[6]+u[54]))^p[11]/((0.5*(u[5]+u[6]+u[54]))^p[11]+(p[23]*p[12])^p[11])*(1-u[8]^p[13]/(u[8]^p[13]+(p[22]*p[14])^p[13])))
        rate[57] = p[21]*(p[27]+p[31]*(0.5*(u[5]+u[6]+u[54]))^p[11]/((0.5*(u[5]+u[6]+u[54]))^p[11]+(p[23]*p[12])^p[11])*(1-u[55]^p[13]/(u[55]^p[13]+(p[22]*p[14])^p[13])))
        # cXR production (scaling factor * (1-silencing))
        rate[3] = p[22]*(1-u[51]^p[5]/(u[51]^p[5]+(p[21]*p[6])^p[5]))
        rate[4] = p[22]*(1-u[49]^p[5]/(u[49]^p[5]+(p[21]*p[6])^p[5]))
        rate[58] = p[22]*(1-u[97]^p[5]/(u[97]^p[5]+(p[21]*p[6])^p[5]))
        # tXA production (scaling factor * (1-silencing))
        rate[5] = p[23]*2/p[33]*(1-u[52]^p[3]/(u[52]^p[3]+(p[21]*p[4])^p[3]))
        rate[6] = p[23]*2/p[33]*(1-u[50]^p[3]/(u[50]^p[3]+(p[21]*p[4])^p[3]))
        rate[59] = p[23]*2/p[33]*(1-u[98]^p[3]/(u[98]^p[3]+(p[21]*p[4])^p[3]))
        # degradation
        rate[7] = 0.1733*u[3]
        rate[8] = 0.1733*u[4]
        rate[60] = 0.1733*u[101]
        rate[9] = 1*u[7]
        rate[10] = 1*u[8]
        rate[61] = 1*u[55]
        rate[11] = 1*u[5]
        rate[12] = 1*u[6]
        rate[62] = 1*u[54]
        rate[13] = 1*u[49]
        rate[14] = 1*u[50]
        rate[63] = 1*u[97]
        rate[15] = 1*u[51]
        rate[16] = 1*u[52]
        rate[64] = 1*u[98]
        ### silencing intermediates on chr1
        rate[17] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[1],0)
        rate[18] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[10],0)
        rate[19] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[11],0)
        rate[20] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[12],0)
        rate[21] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[13],0)
        rate[22] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[14],0)
        rate[23] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[15],0)
        rate[24] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[16],0)
        rate[25] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[17],0)
        rate[26] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[18],0)
        rate[27] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[19],0)
        rate[28] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[20],0)
        rate[29] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[21],0)
        rate[30] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[22],0)
        rate[31] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[23],0)
        rate[32] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[24],0)
        rate[33] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[25],0)
        rate[34] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[26],0)
        rate[35] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[27],0)
        rate[36] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[28],0)
        ### silencing intermediates chr 2
        rate[37] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[2],0)
        rate[38] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[30],0)
        rate[39] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[31],0)
        rate[40] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[32],0)
        rate[41] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[33],0)
        rate[42] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[34],0)
        rate[43] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[35],0)
        rate[44] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[36],0)
        rate[45] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[37],0)
        rate[46] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[38],0)
        rate[47] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[39],0)
        rate[48] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[40],0)
        rate[49] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[41],0)
        rate[50] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[42],0)
        rate[51] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[43],0)
        rate[52] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[44],0)
        rate[53] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[45],0)
        rate[54] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[46],0)
        rate[55] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[47],0)
        rate[56] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[48],0)
        ### silencing intermediates chr 3
        rate[65] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[53],0)
        rate[66] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[59],0)
        rate[67] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[60],0)
        rate[68] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[61],0)
        rate[69] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[62],0)
        rate[70] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[63],0)
        rate[71] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[64],0)
        rate[72] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[65],0)
        rate[73] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[66],0)
        rate[74] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[67],0)
        rate[75] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[68],0)
        rate[76] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[69],0)
        rate[77] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[70],0)
        rate[78] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[71],0)
        rate[79] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[72],0)
        rate[80] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[73],0)
        rate[81] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[74],0)
        rate[82] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[75],0)
        rate[83] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[76],0)
        rate[84] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[77],0)
    elseif p[26] == 4
        # Xist transcription
        rate[1] = p[21]*(p[24]+p[29]*(0.5*(u[5]+u[6]+u[54]+u[57]))^p[11]/((0.5*(u[5]+u[6]+u[54]+u[57]))^p[11]+(p[23]*p[12])^p[11])*(1-u[7]^p[13]/(u[7]^p[13]+(p[22]*p[14])^p[13])))
        rate[2] = p[21]*(p[25]+p[30]*(0.5*(u[5]+u[6]+u[54]+u[57]))^p[11]/((0.5*(u[5]+u[6]+u[54]+u[57]))^p[11]+(p[23]*p[12])^p[11])*(1-u[8]^p[13]/(u[8]^p[13]+(p[22]*p[14])^p[13])))
        rate[57] = p[21]*(p[27]+p[31]*(0.5*(u[5]+u[6]+u[54]+u[57]))^p[11]/((0.5*(u[5]+u[6]+u[54]+u[57]))^p[11]+(p[23]*p[12])^p[11])*(1-u[55]^p[13]/(u[55]^p[13]+(p[22]*p[14])^p[13])))
        rate[85] = p[21]*(p[28]+p[32]*(0.5*(u[5]+u[6]+u[54]+u[57]))^p[11]/((0.5*(u[5]+u[6]+u[54]+u[57]))^p[11]+(p[23]*p[12])^p[11])*(1-u[58]^p[13]/(u[58]^p[13]+(p[22]*p[14])^p[13])))
        # cXR production (scaling factor * (1-silencing))
        rate[3] = p[22]*(1-u[51]^p[5]/(u[51]^p[5]+(p[21]*p[6])^p[5]))
        rate[4] = p[22]*(1-u[49]^p[5]/(u[49]^p[5]+(p[21]*p[6])^p[5]))
        rate[58] = p[22]*(1-u[97]^p[5]/(u[97]^p[5]+(p[21]*p[6])^p[5]))
        rate[86] = p[22]*(1-u[99]^p[5]/(u[99]^p[5]+(p[21]*p[6])^p[5]))
        # tXA production (scaling factor * (1-silencing))
        rate[5] = p[23]*2/p[33]*(1-u[52]^p[3]/(u[52]^p[3]+(p[21]*p[4])^p[3]))
        rate[6] = p[23]*2/p[33]*(1-u[50]^p[3]/(u[50]^p[3]+(p[21]*p[4])^p[3]))
        rate[59] = p[23]*2/p[33]*(1-u[98]^p[3]/(u[98]^p[3]+(p[21]*p[4])^p[3]))
        rate[87] = p[23]*2/p[33]*(1-u[100]^p[3]/(u[100]^p[3]+(p[21]*p[4])^p[3]))
        # degradation
        rate[7] = 0.1733*u[3]
        rate[8] = 0.1733*u[4]
        rate[60] = 0.1733*u[101]
        rate[88] = 0.1733*u[102]
        rate[9] = 1*u[7]
        rate[10] = 1*u[8]
        rate[61] = 1*u[55]
        rate[89] = 1*u[58]
        rate[11] = 1*u[5]
        rate[12] = 1*u[6]
        rate[62] = 1*u[54]
        rate[90] = 1*u[57]
        rate[13] = 1*u[49]
        rate[14] = 1*u[50]
        rate[63] = 1*u[97]
        rate[91] = 1*u[99]
        rate[15] = 1*u[51]
        rate[16] = 1*u[52]
        rate[64] = 1*u[98]
        rate[92] = 1*u[100]
        ### silencing intermediates on chr1
        rate[17] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[1],0)
        rate[18] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[10],0)
        rate[19] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[11],0)
        rate[20] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[12],0)
        rate[21] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[13],0)
        rate[22] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[14],0)
        rate[23] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[15],0)
        rate[24] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[16],0)
        rate[25] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[17],0)
        rate[26] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[18],0)
        rate[27] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[19],0)
        rate[28] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[20],0)
        rate[29] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[21],0)
        rate[30] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[22],0)
        rate[31] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[23],0)
        rate[32] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[24],0)
        rate[33] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[25],0)
        rate[34] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[26],0)
        rate[35] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[27],0)
        rate[36] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[28],0)
        ### silencing intermediates chr 2
        rate[37] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[2],0)
        rate[38] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[30],0)
        rate[39] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[31],0)
        rate[40] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[32],0)
        rate[41] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[33],0)
        rate[42] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[34],0)
        rate[43] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[35],0)
        rate[44] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[36],0)
        rate[45] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[37],0)
        rate[46] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[38],0)
        rate[47] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[39],0)
        rate[48] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[40],0)
        rate[49] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[41],0)
        rate[50] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[42],0)
        rate[51] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[43],0)
        rate[52] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[44],0)
        rate[53] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[45],0)
        rate[54] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[46],0)
        rate[55] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[47],0)
        rate[56] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[48],0)
        ### silencing intermediates chr 3
        rate[65] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[53],0)
        rate[66] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[59],0)
        rate[67] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[60],0)
        rate[68] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[61],0)
        rate[69] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[62],0)
        rate[70] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[63],0)
        rate[71] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[64],0)
        rate[72] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[65],0)
        rate[73] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[66],0)
        rate[74] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[67],0)
        rate[75] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[68],0)
        rate[76] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[69],0)
        rate[77] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[70],0)
        rate[78] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[71],0)
        rate[79] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[72],0)
        rate[80] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[73],0)
        rate[81] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[74],0)
        rate[82] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[75],0)
        rate[83] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[76],0)
        rate[84] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[77],0)
        ### silencing intermediates chr 4
        rate[93] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[56],0)
        rate[94] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[78],0)
        rate[95] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[79],0)
        rate[96] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[80],0)
        rate[97] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[81],0)
        rate[98] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[82],0)
        rate[99] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[83],0)
        rate[100] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[84],0)
        rate[101] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[85],0)
        rate[102] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[86],0)
        rate[103] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[87],0)
        rate[104] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[88],0)
        rate[105] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[89],0)
        rate[106] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[90],0)
        rate[107] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[91],0)
        rate[108] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[92],0)
        rate[109] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[93],0)
        rate[110] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[94],0)
        rate[111] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[95],0)
        rate[112] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[96],0)
    end
	
	# Gillespie
    r = rand(2)
    # determine time point of next reaction
    t_step::Float64 = 1/sum(rate)*log(1/(1-r[1]))
    # decide which reaction will occur
    temp::Float64 = r[2]*sum(rate)
    cum_rate = cumsum(rate)
    # Execute reaction
    if temp<cum_rate[1]
      u[1] += 1
      u[3] += 1
        if p[7]==0
        u[49] += 1
      end
      if p[8]==0
        u[50] += 1
      end
  elseif temp<cum_rate[2]
      u[2] += 1
      u[4] += 1
      if p[7]==0
        u[51] += 1
      end
      if p[8]==0
        u[52] += 1
      end
  elseif temp<cum_rate[3]
        u[8] += 1
  elseif temp<cum_rate[4]
      u[7] += 1
  elseif temp<cum_rate[5]
        u[6] += 1
  elseif temp<cum_rate[6]
        u[5] += 1
  elseif temp<cum_rate[7]
      if maximum([p[7],p[8]])==0
        u[1] -= 1
      end
        u[3] -=1
  elseif temp<cum_rate[8]
      if maximum([p[7],p[8]])==0
        u[2] -= 1
      end
        u[4] -=1
  elseif temp<cum_rate[9]
        u[7] -= 1
  elseif temp<cum_rate[10]
        u[8] -= 1
  elseif temp<cum_rate[11]
        u[5] -= 1
  elseif temp<cum_rate[12]
        u[6] -= 1
  elseif temp<cum_rate[13]
        u[49] -= 1
  elseif temp<cum_rate[14]
        u[50] -= 1
  elseif temp<cum_rate[15]
        u[51] -= 1
  elseif temp<cum_rate[16]
        u[52] -= 1
  elseif temp<cum_rate[17]
      if maximum([p[7],p[8]])>1
        u[10] += 1
        u[1] -= 1
      else
        u[1] -= 1
      end
      p[7]==1 ? u[49] += 1 : 0
      p[8]==1 ? u[50] += 1 : 0
  elseif temp<cum_rate[18]
      if maximum([p[7],p[8]])>2
        u[10] -= 1
        u[11] += 1
      else
        u[10] -= 1
      end
      p[7]==2 ? u[49] += 1 : 0
      p[8]==2 ? u[50] += 1 : 0
  elseif temp<cum_rate[19]
      if maximum([p[7],p[8]])>3
        u[11] -= 1
        u[12] += 1
      else
        u[11] -= 1
      end
      p[7]==3 ? u[49] += 1 : 0
      p[8]==3 ? u[50] += 1 : 0
  elseif temp<cum_rate[20]
      if maximum([p[7],p[8]])>4
        u[12] -= 1
        u[13] += 1
      else
        u[12] -= 1
      end
      p[7]==4 ? u[49] += 1 : 0
      p[8]==4 ? u[50] += 1 : 0
  elseif temp<cum_rate[21]
      if maximum([p[7],p[8]])>5
        u[13] -= 1
        u[14] += 1
      else
        u[13] -= 1
      end
      p[7]==5 ? u[49] += 1 : 0
      p[8]==5 ? u[50] += 1 : 0
  elseif temp<cum_rate[22]
      if maximum([p[7],p[8]])>6
        u[14] -= 1
        u[15] += 1
      else
        u[14] -= 1
      end
      p[7]==6 ? u[49] += 1 : 0
      p[8]==6 ? u[50] += 1 : 0
  elseif temp<cum_rate[23]
      if maximum([p[7],p[8]])>7
        u[15] -= 1
        u[16] += 1
      else
        u[15] -= 1
      end
      p[7]==7 ? u[49] += 1 : 0
      p[8]==7 ? u[50] += 1 : 0
  elseif temp<cum_rate[24]
      if maximum([p[7],p[8]])>8
        u[16] -= 1
        u[17] += 1
      else
        u[16] -= 1
      end
      p[7]==8 ? u[49] += 1 : 0
      p[8]==8 ? u[50] += 1 : 0
  elseif temp<cum_rate[25]
      if maximum([p[7],p[8]])>9
        u[17] -= 1
        u[18] += 1
      else
        u[17] -= 1
      end
      p[7]==9 ? u[49] += 1 : 0
      p[8]==9 ? u[50] += 1 : 0
  elseif temp<cum_rate[26]
      if maximum([p[7],p[8]])>10
        u[18] -= 1
        u[19] += 1
      else
        u[18] -= 1
      end
      p[7]==10 ? u[49] += 1 : 0
      p[8]==10 ? u[50] += 1 : 0
  elseif temp<cum_rate[27]
      if maximum([p[7],p[8]])>11
        u[19] -= 1
        u[20] += 1
      else
        u[19] -= 1
      end
      p[7]==11 ? u[49] += 1 : 0
      p[8]==11 ? u[50] += 1 : 0
  elseif temp<cum_rate[28]
      if maximum([p[7],p[8]])>12
        u[20] -= 1
        u[21] += 1
      else
        u[20] -= 1
      end
      p[7]==12 ? u[49] += 1 : 0
      p[8]==12 ? u[50] += 1 : 0
  elseif temp<cum_rate[29]
      if maximum([p[7],p[8]])>13
        u[21] -= 1
        u[22] += 1
      else
        u[21] -= 1
      end
      p[7]==13 ? u[49] += 1 : 0
      p[8]==13 ? u[50] += 1 : 0
  elseif temp<cum_rate[30]
      if maximum([p[7],p[8]])>14
        u[22] -= 1
        u[23] += 1
      else
        u[22] -= 1
      end
      p[7]==14 ? u[49] += 1 : 0
      p[8]==14 ? u[50] += 1 : 0
  elseif temp<cum_rate[31]
      if maximum([p[7],p[8]])>15
        u[23] -= 1
        u[24] += 1
      else
        u[23] -= 1
      end
      p[7]==15 ? u[49] += 1 : 0
      p[8]==15 ? u[50] += 1 : 0
  elseif temp<cum_rate[32]
      if maximum([p[7],p[8]])>16
        u[24] -= 1
        u[25] += 1
      else
        u[24] -= 1
      end
      p[7]==16 ? u[49] += 1 : 0
      p[8]==16 ? u[50] += 1 : 0
  elseif temp<cum_rate[33]
      if maximum([p[7],p[8]])>17
        u[25] -= 1
        u[26] += 1
      else
        u[25] -= 1
      end
      p[7]==17 ? u[49] += 1 : 0
      p[8]==17 ? u[50] += 1 : 0
  elseif temp<cum_rate[34]
      if maximum([p[7],p[8]])>18
        u[26] -= 1
        u[27] += 1
      else
        u[26] -= 1
      end
      p[7]==18 ? u[49] += 1 : 0
      p[8]==18 ? u[50] += 1 : 0
  elseif temp<cum_rate[35]
      if maximum([p[7],p[8]])>19
        u[27] -= 1
        u[28] += 1
      else
        u[27] -= 1
      end
      p[7]==19 ? u[49] += 1 : 0
      p[8]==19 ? u[50] += 1 : 0
  elseif temp<cum_rate[36]
      u[28] -= 1
      p[7]==20 ? u[49] += 1 : 0
      p[8]==20 ? u[50] += 1 : 0
  elseif temp<cum_rate[37]
      if maximum([p[7],p[8]])>1
        u[30] += 1
        u[2] -= 1
      else
        u[2] -=1
      end
      p[7]==1 ? u[51] += 1 : 0
      p[8]==1 ? u[52] += 1 : 0
  elseif temp<cum_rate[38]
      if maximum([p[7],p[8]])>2
        u[30] -= 1
        u[31] += 1
      else
        u[30] -= 1
      end
      p[7]==2 ? u[51] += 1 : 0
      p[8]==2 ? u[52] += 1 : 0
  elseif temp<cum_rate[39]
      if maximum([p[7],p[8]])>3
        u[31] -= 1
        u[32] += 1
      else
        u[31] -= 1
      end
      p[7]==3 ? u[51] += 1 : 0
      p[8]==3 ? u[52] += 1 : 0
  elseif temp<cum_rate[40]
      if maximum([p[7],p[8]])>4
        u[32] -= 1
        u[33] += 1
      else
        u[32] -= 1
      end
      p[7]==4 ? u[51] += 1 : 0
      p[8]==4 ? u[52] += 1 : 0
  elseif temp<cum_rate[41]
      if maximum([p[7],p[8]])>5
        u[33] -= 1
        u[34] += 1
      else
        u[33] -= 1
      end
      p[7]==5 ? u[51] += 1 : 0
      p[8]==5 ? u[52] += 1 : 0
  elseif temp<cum_rate[42]
      if maximum([p[7],p[8]])>6
        u[34] -= 1
        u[35] += 1
      else
        u[34] -= 1
      end
      p[7]==6 ? u[51] += 1 : 0
      p[8]==6 ? u[52] += 1 : 0
  elseif temp<cum_rate[43]
      if maximum([p[7],p[8]])>7
        u[35] -= 1
        u[36] += 1
      else
        u[35] -= 1
      end
      p[7]==7 ? u[51] += 1 : 0
      p[8]==7 ? u[52] += 1 : 0
  elseif temp<cum_rate[44]
      if maximum([p[7],p[8]])>8
        u[36] -= 1
        u[37] += 1
      else
        u[36] -= 1
      end
      p[7]==8 ? u[51] += 1 : 0
      p[8]==8 ? u[52] += 1 : 0
  elseif temp<cum_rate[45]
      if maximum([p[7],p[8]])>9
        u[37] -= 1
        u[38] += 1
      else
        u[37] -= 1
      end
      p[7]==9 ? u[51] += 1 : 0
      p[8]==9 ? u[52] += 1 : 0
  elseif temp<cum_rate[46]
      if maximum([p[7],p[8]])>10
        u[38] -= 1
        u[39] += 1
      else
        u[38] -= 1
      end
      p[7]==10 ? u[51] += 1 : 0
      p[8]==10 ? u[52] += 1 : 0
  elseif temp<cum_rate[47]
      if maximum([p[7],p[8]])>11
        u[39] -= 1
        u[40] += 1
      else
        u[39] -= 1
      end
      p[7]==11 ? u[51] += 1 : 0
      p[8]==11 ? u[52] += 1 : 0
  elseif temp<cum_rate[48]
      if maximum([p[7],p[8]])>12
        u[40] -= 1
        u[41] += 1
      else
        u[40] -= 1
      end
      p[7]==12 ? u[51] += 1 : 0
      p[8]==12 ? u[52] += 1 : 0
  elseif temp<cum_rate[49]
      if maximum([p[7],p[8]])>13
        u[41] -= 1
        u[42] += 1
      else
        u[41] -= 1
      end
      p[7]==13 ? u[51] += 1 : 0
      p[8]==13 ? u[52] += 1 : 0
  elseif temp<cum_rate[50]
      if maximum([p[7],p[8]])>14
        u[42] -= 1
        u[43] += 1
      else
        u[42] -= 1
      end
      p[7]==14 ? u[51] += 1 : 0
      p[8]==14 ? u[52] += 1 : 0
  elseif temp<cum_rate[51]
      if maximum([p[7],p[8]])>15
        u[43] -= 1
        u[44] += 1
      else
        u[43] -= 1
      end
      p[7]==15 ? u[51] += 1 : 0
      p[8]==15 ? u[52] += 1 : 0
  elseif temp<cum_rate[52]
      if maximum([p[7],p[8]])>16
        u[44] -= 1
        u[45] += 1
      else
        u[44] -= 1
      end
      p[7]==16 ? u[51] += 1 : 0
      p[8]==16 ? u[52] += 1 : 0
  elseif temp<cum_rate[53]
      if maximum([p[7],p[8]])>17
        u[45] -= 1
        u[46] += 1
      else
        u[45] -= 1
      end
      p[7]==17 ? u[51] += 1 : 0
      p[8]==17 ? u[52] += 1 : 0
  elseif temp<cum_rate[54]
      if maximum([p[7],p[8]])>18
        u[46] -= 1
        u[47] += 1
      else
        u[46] -= 1
      end
      p[7]==18 ? u[51] += 1 : 0
      p[8]==18 ? u[52] += 1 : 0
  elseif temp<cum_rate[55]
      if maximum([p[7],p[8]])>19
        u[47] -= 1
        u[48] += 1
      else
        u[47] -= 1
      end
      p[7]==19 ? u[51] += 1 : 0
      p[8]==19 ? u[52] += 1 : 0
  elseif temp<cum_rate[56]
      u[48] -= 1
      p[7]==20 ? u[51] += 1 : 0
      p[8]==20 ? u[52] += 1 : 0
  elseif temp<cum_rate[57]
      u[53] += 1
      u[101] += 1
      if p[7]==0
        u[97] += 1
      end
      if p[8]==0
        u[98] += 1
      end
  elseif temp<cum_rate[58]
      u[55] += 1
  elseif temp<cum_rate[59]
      u[54] += 1
  elseif temp<cum_rate[60]
      if maximum([p[7],p[8]])==0
        u[53] -= 1
      end
        u[101] -=1
  elseif temp<cum_rate[61]
      u[55] -= 1
  elseif temp<cum_rate[62]
      u[54] -= 1
  elseif temp<cum_rate[63]
      u[97] -= 1
  elseif temp<cum_rate[64]
      u[98] -= 1
  elseif temp<cum_rate[65]
      if maximum([p[7],p[8]])>1
        u[59] += 1
        u[53] -= 1
      else
        u[53] -=1
      end
      p[7]==1 ? u[97] += 1 : 0
      p[8]==1 ? u[98] += 1 : 0
  elseif temp<cum_rate[66]
      if maximum([p[7],p[8]])>2
        u[59] -= 1
        u[60] += 1
      else
        u[59] -= 1
      end
      p[7]==2 ? u[97] += 1 : 0
      p[8]==2 ? u[98] += 1 : 0
  elseif temp<cum_rate[67]
      if maximum([p[7],p[8]])>3
        u[60] -= 1
        u[61] += 1
      else
        u[60] -= 1
      end
      p[7]==3 ? u[97] += 1 : 0
      p[8]==3 ? u[98] += 1 : 0
  elseif temp<cum_rate[68]
      if maximum([p[7],p[8]])>4
        u[61] -= 1
        u[62] += 1
      else
        u[61] -= 1
      end
      p[7]==4 ? u[97] += 1 : 0
      p[8]==4 ? u[98] += 1 : 0
  elseif temp<cum_rate[69]
      if maximum([p[7],p[8]])>5
        u[62] -= 1
        u[63] += 1
      else
        u[62] -= 1
      end
      p[7]==5 ? u[97] += 1 : 0
      p[8]==5 ? u[98] += 1 : 0
  elseif temp<cum_rate[70]
      if maximum([p[7],p[8]])>6
        u[63] -= 1
        u[64] += 1
      else
        u[63] -= 1
      end
      p[7]==6 ? u[97] += 1 : 0
      p[8]==6 ? u[98] += 1 : 0
  elseif temp<cum_rate[71]
      if maximum([p[7],p[8]])>7
        u[64] -= 1
        u[65] += 1
      else
        u[64] -= 1
      end
      p[7]==7 ? u[97] += 1 : 0
      p[8]==7 ? u[98] += 1 : 0
  elseif temp<cum_rate[72]
      if maximum([p[7],p[8]])>8
        u[65] -= 1
        u[66] += 1
      else
        u[65] -= 1
      end
      p[7]==8 ? u[97] += 1 : 0
      p[8]==8 ? u[98] += 1 : 0
  elseif temp<cum_rate[73]
      if maximum([p[7],p[8]])>9
        u[66] -= 1
        u[67] += 1
      else
        u[66] -= 1
      end
      p[7]==9 ? u[97] += 1 : 0
      p[8]==9 ? u[98] += 1 : 0
  elseif temp<cum_rate[74]
      if maximum([p[7],p[8]])>10
        u[67] -= 1
        u[68] += 1
      else
        u[67] -= 1
      end
      p[7]==10 ? u[97] += 1 : 0
      p[8]==10 ? u[98] += 1 : 0
  elseif temp<cum_rate[75]
      if maximum([p[7],p[8]])>11
        u[68] -= 1
        u[69] += 1
      else
        u[68] -= 1
      end
      p[7]==11 ? u[97] += 1 : 0
      p[8]==11 ? u[98] += 1 : 0
  elseif temp<cum_rate[76]
      if maximum([p[7],p[8]])>12
        u[69] -= 1
        u[70] += 1
      else
        u[69] -= 1
      end
      p[7]==12 ? u[97] += 1 : 0
      p[8]==12 ? u[98] += 1 : 0
  elseif temp<cum_rate[77]
      if maximum([p[7],p[8]])>13
        u[70] -= 1
        u[71] += 1
      else
        u[70] -= 1
      end
      p[7]==13 ? u[97] += 1 : 0
      p[8]==13 ? u[98] += 1 : 0
  elseif temp<cum_rate[78]
      if maximum([p[7],p[8]])>14
        u[71] -= 1
        u[72] += 1
      else
        u[71] -= 1
      end
      p[7]==14 ? u[97] += 1 : 0
      p[8]==14 ? u[98] += 1 : 0
  elseif temp<cum_rate[79]
      if maximum([p[7],p[8]])>15
        u[72] -= 1
        u[73] += 1
      else
        u[72] -= 1
      end
      p[7]==15 ? u[97] += 1 : 0
      p[8]==15 ? u[98] += 1 : 0
  elseif temp<cum_rate[80]
      if maximum([p[7],p[8]])>16
        u[73] -= 1
        u[74] += 1
      else
        u[73] -= 1
      end
      p[7]==16 ? u[97] += 1 : 0
      p[8]==16 ? u[98] += 1 : 0
  elseif temp<cum_rate[81]
      if maximum([p[7],p[8]])>17
        u[74] -= 1
        u[75] += 1
      else
        u[74] -= 1
      end
      p[7]==17 ? u[97] += 1 : 0
      p[8]==17 ? u[98] += 1 : 0
  elseif temp<cum_rate[82]
      if maximum([p[7],p[8]])>18
        u[75] -= 1
        u[76] += 1
      else
        u[75] -= 1
      end
      p[7]==18 ? u[97] += 1 : 0
      p[8]==18 ? u[98] += 1 : 0
  elseif temp<cum_rate[83]
      if maximum([p[7],p[8]])>19
        u[76] -= 1
        u[77] += 1
      else
        u[76] -= 1
      end
      p[7]==19 ? u[97] += 1 : 0
      p[8]==19 ? u[98] += 1 : 0
  elseif temp<cum_rate[84]
      u[77] -= 1
      p[7]==20 ? u[97] += 1 : 0
      p[8]==20 ? u[98] += 1 : 0
  elseif temp<cum_rate[85]
      u[56] += 1
      u[102] += 1
      if p[7]==0
        u[99] += 1
      end
      if p[8]==0
        u[100] += 1
      end
  elseif temp<cum_rate[86]
      u[58] += 1
  elseif temp<cum_rate[87]
      u[57] += 1
  elseif temp<cum_rate[88]
      if maximum([p[7],p[8]])==0
        u[56] -= 1
      end
        u[102] -=1
  elseif temp<cum_rate[89]
      u[58] -= 1
  elseif temp<cum_rate[90]
      u[57] -= 1
  elseif temp<cum_rate[91]
      u[99] -= 1
  elseif temp<cum_rate[92]
      u[100] -= 1
  elseif temp<cum_rate[93]
      if maximum([p[7],p[8]])>1
        u[78] += 1
        u[56] -= 1
      else
        u[56] -=1
      end
      p[7]==1 ? u[99] += 1 : 0
      p[8]==1 ? u[100] += 1 : 0
  elseif temp<cum_rate[94]
      if maximum([p[7],p[8]])>2
        u[78] -= 1
        u[79] += 1
      else
        u[78] -= 1
      end
      p[7]==2 ? u[99] += 1 : 0
      p[8]==2 ? u[100] += 1 : 0
  elseif temp<cum_rate[95]
      if maximum([p[7],p[8]])>3
        u[79] -= 1
        u[80] += 1
      else
        u[79] -= 1
      end
      p[7]==3 ? u[99] += 1 : 0
      p[8]==3 ? u[100] += 1 : 0
  elseif temp<cum_rate[96]
      if maximum([p[7],p[8]])>4
        u[80] -= 1
        u[81] += 1
      else
        u[80] -= 1
      end
      p[7]==4 ? u[99] += 1 : 0
      p[8]==4 ? u[100] += 1 : 0
  elseif temp<cum_rate[97]
      if maximum([p[7],p[8]])>5
        u[81] -= 1
        u[82] += 1
      else
        u[81] -= 1
      end
      p[7]==5 ? u[99] += 1 : 0
      p[8]==5 ? u[100] += 1 : 0
  elseif temp<cum_rate[98]
      if maximum([p[7],p[8]])>6
        u[82] -= 1
        u[83] += 1
      else
        u[82] -= 1
      end
      p[7]==6 ? u[99] += 1 : 0
      p[8]==6 ? u[100] += 1 : 0
  elseif temp<cum_rate[99]
      if maximum([p[7],p[8]])>7
        u[83] -= 1
        u[84] += 1
      else
        u[83] -= 1
      end
      p[7]==7 ? u[99] += 1 : 0
      p[8]==7 ? u[100] += 1 : 0
  elseif temp<cum_rate[100]
      if maximum([p[7],p[8]])>8
        u[84] -= 1
        u[85] += 1
      else
        u[84] -= 1
      end
      p[7]==8 ? u[99] += 1 : 0
      p[8]==8 ? u[100] += 1 : 0
  elseif temp<cum_rate[101]
      if maximum([p[7],p[8]])>9
        u[85] -= 1
        u[86] += 1
      else
        u[85] -= 1
      end
      p[7]==9 ? u[99] += 1 : 0
      p[8]==9 ? u[100] += 1 : 0
  elseif temp<cum_rate[102]
      if maximum([p[7],p[8]])>10
        u[86] -= 1
        u[87] += 1
      else
        u[86] -= 1
      end
      p[7]==10 ? u[99] += 1 : 0
      p[8]==10 ? u[100] += 1 : 0
  elseif temp<cum_rate[103]
      if maximum([p[7],p[8]])>11
        u[87] -= 1
        u[88] += 1
      else
        u[87] -= 1
      end
      p[7]==11 ? u[99] += 1 : 0
      p[8]==11 ? u[100] += 1 : 0
  elseif temp<cum_rate[104]
      if maximum([p[7],p[8]])>12
        u[88] -= 1
        u[89] += 1
      else
        u[88] -= 1
      end
      p[7]==12 ? u[99] += 1 : 0
      p[8]==12 ? u[100] += 1 : 0
  elseif temp<cum_rate[105]
      if maximum([p[7],p[8]])>13
        u[89] -= 1
        u[90] += 1
      else
        u[89] -= 1
      end
      p[7]==13 ? u[99] += 1 : 0
      p[8]==13 ? u[100] += 1 : 0
  elseif temp<cum_rate[106]
      if maximum([p[7],p[8]])>14
        u[90] -= 1
        u[91] += 1
      else
        u[90] -= 1
      end
      p[7]==14 ? u[99] += 1 : 0
      p[8]==14 ? u[100] += 1 : 0
  elseif temp<cum_rate[107]
      if maximum([p[7],p[8]])>15
        u[91] -= 1
        u[92] += 1
      else
        u[91] -= 1
      end
      p[7]==15 ? u[99] += 1 : 0
      p[8]==15 ? u[100] += 1 : 0
  elseif temp<cum_rate[108]
      if maximum([p[7],p[8]])>16
        u[92] -= 1
        u[93] += 1
      else
        u[92] -= 1
      end
      p[7]==16 ? u[99] += 1 : 0
      p[8]==16 ? u[100] += 1 : 0
  elseif temp<cum_rate[109]
      if maximum([p[7],p[8]])>17
        u[93] -= 1
        u[94] += 1
      else
        u[93] -= 1
      end
      p[7]==17 ? u[99] += 1 : 0
      p[8]==17 ? u[100] += 1 : 0
  elseif temp<cum_rate[110]
      if maximum([p[7],p[8]])>18
        u[94] -= 1
        u[95] += 1
      else
        u[94] -= 1
      end
      p[7]==18 ? u[99] += 1 : 0
      p[8]==18 ? u[100] += 1 : 0
  elseif temp<cum_rate[111]
      if maximum([p[7],p[8]])>19
        u[95] -= 1
        u[96] += 1
      else
        u[95] -= 1
      end
      p[7]==19 ? u[99] += 1 : 0
      p[8]==19 ? u[100] += 1 : 0
  elseif temp<cum_rate[112]
      u[96] -= 1
      p[7]==20 ? u[99] += 1 : 0
      p[8]==20 ? u[100] += 1 : 0
  end
    t = t+t_step
    push!(out,copy(u))
    push!(t_out,copy(t))
  end
  return out, t_out
end

### Dampening with cXR upregulation hypothesis
function sim_delay_damp_cXR_upreg(p::Array{Float64,1},tspan::Tuple{Float64,Float64},u0::Array{Int64,1})
  t::Float64=tspan[1]
  i::Int64=0
  out::Array{Array{Int64,1}} = []
  rate_out::Array{Array{Float64,1}} = []
  t_out::Array{Float64,1} = []
  push!(out, u0)
  push!(t_out,t)
  u::Array{Int64,1}=copy(u0)
  #############################################################
  rate::Array{Float64,1}=zeros(56)
  while (t<tspan[2])
    i+=1
    p[24] = 0
    p[29] = 1
    p[25] = 0
    p[30] = 1
    # Xist transcription
    rate[1] = p[21]*(p[24]+p[29]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[7]^p[13]/(u[7]^p[13]+(p[22]*p[14])^p[13])))
    rate[2] = p[21]*(p[25]+p[30]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[8]^p[13]/(u[8]^p[13]+(p[22]*p[14])^p[13])))
    # cXR production
    if t>100
        rate[3] = p[22]*(1-u[51]^p[5]/(u[51]^p[5]+(p[21]*p[6])^p[5]))
        rate[4] = p[22]*(1-u[49]^p[5]/(u[49]^p[5]+(p[21]*p[6])^p[5]))
    end
    # tXA production
    if t<=100
        rate[5] = p[23]
        rate[6] = p[23]
    else
        rate[5] = p[23]*(1-u[52]^p[3]/(u[52]^p[3]+(p[21]*p[4])^p[3]))
        rate[6] = p[23]*(1-u[50]^p[3]/(u[50]^p[3]+(p[21]*p[4])^p[3]))
    end
    # degradation
    rate[7] = 0.1733*u[3]
    rate[8] = 0.1733*u[4]
    rate[9] = 1*u[7]
    rate[10] = 1*u[8]
    rate[11] = 1*u[5]
    rate[12] = 1*u[6]
    rate[13] = 1*u[49]
    rate[14] = 1*u[50]
    rate[15] = 1*u[51]
    rate[16] = 1*u[52]
    ### silencing intermediates on chr1
    if t>100
        rate[17] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[1],0)
        rate[18] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[10],0)
        rate[19] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[11],0)
        rate[20] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[12],0)
        rate[21] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[13],0)
        rate[22] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[14],0)
        rate[23] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[15],0)
        rate[24] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[16],0)
        rate[25] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[17],0)
        rate[26] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[18],0)
        rate[27] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[19],0)
        rate[28] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[20],0)
        rate[29] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[21],0)
        rate[30] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[22],0)
        rate[31] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[23],0)
        rate[32] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[24],0)
        rate[33] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[25],0)
        rate[34] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[26],0)
        rate[35] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[27],0)
        rate[36] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[28],0)
    ### silencing intermediates chr 2
        rate[37] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[2],0)
        rate[38] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[30],0)
        rate[39] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[31],0)
        rate[40] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[32],0)
        rate[41] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[33],0)
        rate[42] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[34],0)
        rate[43] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[35],0)
        rate[44] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[36],0)
        rate[45] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[37],0)
        rate[46] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[38],0)
        rate[47] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[39],0)
        rate[48] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[40],0)
        rate[49] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[41],0)
        rate[50] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[42],0)
        rate[51] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[43],0)
        rate[52] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[44],0)
        rate[53] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[45],0)
        rate[54] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[46],0)
        rate[55] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[47],0)
        rate[56] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[48],0)
    end
    
    # Gillespie
    r = rand(2)
    # determine time point of next reaction
    t_step::Float64 = 1/sum(rate)*log(1/(1-r[1]))
    # decide which reaction will occur
    temp::Float64 = r[2]*sum(rate)
    cum_rate = cumsum(rate)
    # Execute reaction
    if temp<cum_rate[1]
      u[1] += 1
      u[3] += 1
        if p[7]==0
        u[49] += 1
      end
      if p[8]==0
        u[50] += 1
      end
  elseif temp<cum_rate[2]
      u[2] += 1
      u[4] += 1
      if p[7]==0
        u[51] += 1
      end
      if p[8]==0
        u[52] += 1
      end
  elseif temp<cum_rate[3]
        u[8] += 1
  elseif temp<cum_rate[4]
      u[7] += 1
  elseif temp<cum_rate[5]
        u[6] += 1
  elseif temp<cum_rate[6]
        u[5] += 1
  elseif temp<cum_rate[7]
      if maximum([p[7],p[8]])==0
        u[1] -= 1
      end
        u[3] -=1
  elseif temp<cum_rate[8]
      if maximum([p[7],p[8]])==0
        u[2] -= 1
      end
        u[4] -=1
  elseif temp<cum_rate[9]
        u[7] -= 1
  elseif temp<cum_rate[10]
        u[8] -= 1
  elseif temp<cum_rate[11]
        u[5] -= 1
  elseif temp<cum_rate[12]
        u[6] -= 1
  elseif temp<cum_rate[13]
        u[49] -= 1
  elseif temp<cum_rate[14]
        u[50] -= 1
  elseif temp<cum_rate[15]
        u[51] -= 1
  elseif temp<cum_rate[16]
        u[52] -= 1
  elseif temp<cum_rate[17]
      if maximum([p[7],p[8]])>1
        u[10] += 1
        u[1] -= 1
      else
        u[1] -= 1
      end
      p[7]==1 ? u[49] += 1 : 0
      p[8]==1 ? u[50] += 1 : 0
  elseif temp<cum_rate[18]
      if maximum([p[7],p[8]])>2
        u[10] -= 1
        u[11] += 1
      else
        u[10] -= 1
      end
      p[7]==2 ? u[49] += 1 : 0
      p[8]==2 ? u[50] += 1 : 0
  elseif temp<cum_rate[19]
      if maximum([p[7],p[8]])>3
        u[11] -= 1
        u[12] += 1
      else
        u[11] -= 1
      end
      p[7]==3 ? u[49] += 1 : 0
      p[8]==3 ? u[50] += 1 : 0
  elseif temp<cum_rate[20]
      if maximum([p[7],p[8]])>4
        u[12] -= 1
        u[13] += 1
      else
        u[12] -= 1
      end
      p[7]==4 ? u[49] += 1 : 0
      p[8]==4 ? u[50] += 1 : 0
  elseif temp<cum_rate[21]
      if maximum([p[7],p[8]])>5
        u[13] -= 1
        u[14] += 1
      else
        u[13] -= 1
      end
      p[7]==5 ? u[49] += 1 : 0
      p[8]==5 ? u[50] += 1 : 0
  elseif temp<cum_rate[22]
      if maximum([p[7],p[8]])>6
        u[14] -= 1
        u[15] += 1
      else
        u[14] -= 1
      end
      p[7]==6 ? u[49] += 1 : 0
      p[8]==6 ? u[50] += 1 : 0
  elseif temp<cum_rate[23]
      if maximum([p[7],p[8]])>7
        u[15] -= 1
        u[16] += 1
      else
        u[15] -= 1
      end
      p[7]==7 ? u[49] += 1 : 0
      p[8]==7 ? u[50] += 1 : 0
  elseif temp<cum_rate[24]
      if maximum([p[7],p[8]])>8
        u[16] -= 1
        u[17] += 1
      else
        u[16] -= 1
      end
      p[7]==8 ? u[49] += 1 : 0
      p[8]==8 ? u[50] += 1 : 0
  elseif temp<cum_rate[25]
      if maximum([p[7],p[8]])>9
        u[17] -= 1
        u[18] += 1
      else
        u[17] -= 1
      end
      p[7]==9 ? u[49] += 1 : 0
      p[8]==9 ? u[50] += 1 : 0
  elseif temp<cum_rate[26]
      if maximum([p[7],p[8]])>10
        u[18] -= 1
        u[19] += 1
      else
        u[18] -= 1
      end
      p[7]==10 ? u[49] += 1 : 0
      p[8]==10 ? u[50] += 1 : 0
  elseif temp<cum_rate[27]
      if maximum([p[7],p[8]])>11
        u[19] -= 1
        u[20] += 1
      else
        u[19] -= 1
      end
      p[7]==11 ? u[49] += 1 : 0
      p[8]==11 ? u[50] += 1 : 0
  elseif temp<cum_rate[28]
      if maximum([p[7],p[8]])>12
        u[20] -= 1
        u[21] += 1
      else
        u[20] -= 1
      end
      p[7]==12 ? u[49] += 1 : 0
      p[8]==12 ? u[50] += 1 : 0
  elseif temp<cum_rate[29]
      if maximum([p[7],p[8]])>13
        u[21] -= 1
        u[22] += 1
      else
        u[21] -= 1
      end
      p[7]==13 ? u[49] += 1 : 0
      p[8]==13 ? u[50] += 1 : 0
  elseif temp<cum_rate[30]
      if maximum([p[7],p[8]])>14
        u[22] -= 1
        u[23] += 1
      else
        u[22] -= 1
      end
      p[7]==14 ? u[49] += 1 : 0
      p[8]==14 ? u[50] += 1 : 0
  elseif temp<cum_rate[31]
      if maximum([p[7],p[8]])>15
        u[23] -= 1
        u[24] += 1
      else
        u[23] -= 1
      end
      p[7]==15 ? u[49] += 1 : 0
      p[8]==15 ? u[50] += 1 : 0
  elseif temp<cum_rate[32]
      if maximum([p[7],p[8]])>16
        u[24] -= 1
        u[25] += 1
      else
        u[24] -= 1
      end
      p[7]==16 ? u[49] += 1 : 0
      p[8]==16 ? u[50] += 1 : 0
  elseif temp<cum_rate[33]
      if maximum([p[7],p[8]])>17
        u[25] -= 1
        u[26] += 1
      else
        u[25] -= 1
      end
      p[7]==17 ? u[49] += 1 : 0
      p[8]==17 ? u[50] += 1 : 0
  elseif temp<cum_rate[34]
      if maximum([p[7],p[8]])>18
        u[26] -= 1
        u[27] += 1
      else
        u[26] -= 1
      end
      p[7]==18 ? u[49] += 1 : 0
      p[8]==18 ? u[50] += 1 : 0
  elseif temp<cum_rate[35]
      if maximum([p[7],p[8]])>19
        u[27] -= 1
        u[28] += 1
      else
        u[27] -= 1
      end
      p[7]==19 ? u[49] += 1 : 0
      p[8]==19 ? u[50] += 1 : 0
  elseif temp<cum_rate[36]
      u[28] -= 1
      p[7]==20 ? u[49] += 1 : 0
      p[8]==20 ? u[50] += 1 : 0
  elseif temp<cum_rate[37]
      if maximum([p[7],p[8]])>1
        u[30] += 1
        u[2] -= 1
      else
        u[2] -=1
      end
      p[7]==1 ? u[51] += 1 : 0
      p[8]==1 ? u[52] += 1 : 0
  elseif temp<cum_rate[38]
      if maximum([p[7],p[8]])>2
        u[30] -= 1
        u[31] += 1
      else
        u[30] -= 1
      end
      p[7]==2 ? u[51] += 1 : 0
      p[8]==2 ? u[52] += 1 : 0
  elseif temp<cum_rate[39]
      if maximum([p[7],p[8]])>3
        u[31] -= 1
        u[32] += 1
      else
        u[31] -= 1
      end
      p[7]==3 ? u[51] += 1 : 0
      p[8]==3 ? u[52] += 1 : 0
  elseif temp<cum_rate[40]
      if maximum([p[7],p[8]])>4
        u[32] -= 1
        u[33] += 1
      else
        u[32] -= 1
      end
      p[7]==4 ? u[51] += 1 : 0
      p[8]==4 ? u[52] += 1 : 0
  elseif temp<cum_rate[41]
      if maximum([p[7],p[8]])>5
        u[33] -= 1
        u[34] += 1
      else
        u[33] -= 1
      end
      p[7]==5 ? u[51] += 1 : 0
      p[8]==5 ? u[52] += 1 : 0
  elseif temp<cum_rate[42]
      if maximum([p[7],p[8]])>6
        u[34] -= 1
        u[35] += 1
      else
        u[34] -= 1
      end
      p[7]==6 ? u[51] += 1 : 0
      p[8]==6 ? u[52] += 1 : 0
  elseif temp<cum_rate[43]
      if maximum([p[7],p[8]])>7
        u[35] -= 1
        u[36] += 1
      else
        u[35] -= 1
      end
      p[7]==7 ? u[51] += 1 : 0
      p[8]==7 ? u[52] += 1 : 0
  elseif temp<cum_rate[44]
      if maximum([p[7],p[8]])>8
        u[36] -= 1
        u[37] += 1
      else
        u[36] -= 1
      end
      p[7]==8 ? u[51] += 1 : 0
      p[8]==8 ? u[52] += 1 : 0
  elseif temp<cum_rate[45]
      if maximum([p[7],p[8]])>9
        u[37] -= 1
        u[38] += 1
      else
        u[37] -= 1
      end
      p[7]==9 ? u[51] += 1 : 0
      p[8]==9 ? u[52] += 1 : 0
  elseif temp<cum_rate[46]
      if maximum([p[7],p[8]])>10
        u[38] -= 1
        u[39] += 1
      else
        u[38] -= 1
      end
      p[7]==10 ? u[51] += 1 : 0
      p[8]==10 ? u[52] += 1 : 0
  elseif temp<cum_rate[47]
      if maximum([p[7],p[8]])>11
        u[39] -= 1
        u[40] += 1
      else
        u[39] -= 1
      end
      p[7]==11 ? u[51] += 1 : 0
      p[8]==11 ? u[52] += 1 : 0
  elseif temp<cum_rate[48]
      if maximum([p[7],p[8]])>12
        u[40] -= 1
        u[41] += 1
      else
        u[40] -= 1
      end
      p[7]==12 ? u[51] += 1 : 0
      p[8]==12 ? u[52] += 1 : 0
  elseif temp<cum_rate[49]
      if maximum([p[7],p[8]])>13
        u[41] -= 1
        u[42] += 1
      else
        u[41] -= 1
      end
      p[7]==13 ? u[51] += 1 : 0
      p[8]==13 ? u[52] += 1 : 0
  elseif temp<cum_rate[50]
      if maximum([p[7],p[8]])>14
        u[42] -= 1
        u[43] += 1
      else
        u[42] -= 1
      end
      p[7]==14 ? u[51] += 1 : 0
      p[8]==14 ? u[52] += 1 : 0
  elseif temp<cum_rate[51]
      if maximum([p[7],p[8]])>15
        u[43] -= 1
        u[44] += 1
      else
        u[43] -= 1
      end
      p[7]==15 ? u[51] += 1 : 0
      p[8]==15 ? u[52] += 1 : 0
  elseif temp<cum_rate[52]
      if maximum([p[7],p[8]])>16
        u[44] -= 1
        u[45] += 1
      else
        u[44] -= 1
      end
      p[7]==16 ? u[51] += 1 : 0
      p[8]==16 ? u[52] += 1 : 0
  elseif temp<cum_rate[53]
      if maximum([p[7],p[8]])>17
        u[45] -= 1
        u[46] += 1
      else
        u[45] -= 1
      end
      p[7]==17 ? u[51] += 1 : 0
      p[8]==17 ? u[52] += 1 : 0
  elseif temp<cum_rate[54]
      if maximum([p[7],p[8]])>18
        u[46] -= 1
        u[47] += 1
      else
        u[46] -= 1
      end
      p[7]==18 ? u[51] += 1 : 0
      p[8]==18 ? u[52] += 1 : 0
  elseif temp<cum_rate[55]
      if maximum([p[7],p[8]])>19
        u[47] -= 1
        u[48] += 1
      else
        u[47] -= 1
      end
      p[7]==19 ? u[51] += 1 : 0
      p[8]==19 ? u[52] += 1 : 0
  elseif temp<cum_rate[56]
      u[48] -= 1
      p[7]==20 ? u[51] += 1 : 0
      p[8]==20 ? u[52] += 1 : 0
  end
    t = t+t_step
    push!(out,copy(u))
    push!(t_out,copy(t))
  end
  return out, t_out
end

### Dampening with cXR dampening hypothesis
function sim_delay_damp_cXR_damp(p::Array{Float64,1},tspan::Tuple{Float64,Float64},u0::Array{Int64,1})
  t::Float64=tspan[1]
  i::Int64=0
  out::Array{Array{Int64,1}} = []
  t_out::Array{Float64,1} = []
  push!(out, u0)
  push!(t_out,t)
  u::Array{Int64,1}=copy(u0)
  #############################################################
  rate::Array{Float64,1}=zeros(56)
  while (t<tspan[2])
    i+=1
    p[24] = 0
    p[29] = 1
    p[25] = 0
    p[30] = 1
    if t>100
        p[1] = 1
    end
    # Xist transcription
    rate[1] = p[21]*(p[24]+p[29]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[7]^p[13]/(u[7]^p[13]+(p[22]*p[14])^p[13])))
    rate[2] = p[21]*(p[25]+p[30]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[8]^p[13]/(u[8]^p[13]+(p[22]*p[14])^p[13])))
    # cXR production
    if t<100
        rate[3] = p[22]*(1-p[1]*u[2]^p[5]/(u[2]^p[5]+(p[21]*p[6])^p[5]))
        rate[4] = p[22]*(1-p[1]*u[1]^p[5]/(u[1]^p[5]+(p[21]*p[6])^p[5]))
    else
        rate[3] = p[22]*(1-p[1]*u[51]^p[5]/(u[51]^p[5]+(p[21]*p[6])^p[5]))
        rate[4] = p[22]*(1-p[1]*u[49]^p[5]/(u[49]^p[5]+(p[21]*p[6])^p[5]))
    end
    # tXA production
    if t<100
        rate[5] = p[23]
        rate[6] = p[23]
    else
        rate[5] = p[23]*(1-u[52]^p[3]/(u[52]^p[3]+(p[21]*p[4])^p[3]))
        rate[6] = p[23]*(1-u[50]^p[3]/(u[50]^p[3]+(p[21]*p[4])^p[3]))
    end
    # degradation
    rate[7] = 0.1733*u[3]
    rate[8] = 0.1733*u[4]
    rate[9] = 1*u[7]
    rate[10] = 1*u[8]
    rate[11] = 1*u[5]
    rate[12] = 1*u[6]
    rate[13] = 1*u[49]
    rate[14] = 1*u[50]
    rate[15] = 1*u[51]
    rate[16] = 1*u[52]
    ### silencing intermediates on chr1
    if t>100
    rate[17] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[1],0)
    rate[18] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[10],0)
    rate[19] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[11],0)
    rate[20] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[12],0)
    rate[21] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[13],0)
    rate[22] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[14],0)
    rate[23] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[15],0)
    rate[24] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[16],0)
    rate[25] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[17],0)
    rate[26] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[18],0)
    rate[27] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[19],0)
    rate[28] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[20],0)
    rate[29] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[21],0)
    rate[30] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[22],0)
    rate[31] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[23],0)
    rate[32] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[24],0)
    rate[33] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[25],0)
    rate[34] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[26],0)
    rate[35] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[27],0)
    rate[36] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[28],0)
    ### silencing intermediates chr 2
    rate[37] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[2],0)
    rate[38] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[30],0)
    rate[39] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[31],0)
    rate[40] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[32],0)
    rate[41] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[33],0)
    rate[42] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[34],0)
    rate[43] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[35],0)
    rate[44] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[36],0)
    rate[45] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[37],0)
    rate[46] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[38],0)
    rate[47] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[39],0)
    rate[48] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[40],0)
    rate[49] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[41],0)
    rate[50] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[42],0)
    rate[51] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[43],0)
    rate[52] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[44],0)
    rate[53] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[45],0)
    rate[54] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[46],0)
    rate[55] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[47],0)
    rate[56] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[48],0)
    end
    
    # Gillespie
    r = rand(2)
    # determine time point of next reaction
    t_step::Float64 = 1/sum(rate)*log(1/(1-r[1]))
    # decide which reaction will occur
    temp::Float64 = r[2]*sum(rate)
    cum_rate = cumsum(rate)
    # Execute reaction
    if temp<cum_rate[1]
      u[1] += 1
      u[3] += 1
        if p[7]==0
        u[49] += 1
      end
      if p[8]==0
        u[50] += 1
      end
  elseif temp<cum_rate[2]
      u[2] += 1
      u[4] += 1
      if p[7]==0
        u[51] += 1
      end
      if p[8]==0
        u[52] += 1
      end
  elseif temp<cum_rate[3]
        u[8] += 1
  elseif temp<cum_rate[4]
      u[7] += 1
  elseif temp<cum_rate[5]
        u[6] += 1
  elseif temp<cum_rate[6]
        u[5] += 1
  elseif temp<cum_rate[7]
      if maximum([p[7],p[8]])==0
        u[1] -= 1
      end
        u[3] -=1
  elseif temp<cum_rate[8]
      if maximum([p[7],p[8]])==0
        u[2] -= 1
      end
        u[4] -=1
  elseif temp<cum_rate[9]
        u[7] -= 1
  elseif temp<cum_rate[10]
        u[8] -= 1
  elseif temp<cum_rate[11]
        u[5] -= 1
  elseif temp<cum_rate[12]
        u[6] -= 1
  elseif temp<cum_rate[13]
        u[49] -= 1
  elseif temp<cum_rate[14]
        u[50] -= 1
  elseif temp<cum_rate[15]
        u[51] -= 1
  elseif temp<cum_rate[16]
        u[52] -= 1
  elseif temp<cum_rate[17]
      if maximum([p[7],p[8]])>1
        u[10] += 1
        u[1] -= 1
      else
        u[1] -= 1
      end
      p[7]==1 ? u[49] += 1 : 0
      p[8]==1 ? u[50] += 1 : 0
  elseif temp<cum_rate[18]
      if maximum([p[7],p[8]])>2
        u[10] -= 1
        u[11] += 1
      else
        u[10] -= 1
      end
      p[7]==2 ? u[49] += 1 : 0
      p[8]==2 ? u[50] += 1 : 0
  elseif temp<cum_rate[19]
      if maximum([p[7],p[8]])>3
        u[11] -= 1
        u[12] += 1
      else
        u[11] -= 1
      end
      p[7]==3 ? u[49] += 1 : 0
      p[8]==3 ? u[50] += 1 : 0
  elseif temp<cum_rate[20]
      if maximum([p[7],p[8]])>4
        u[12] -= 1
        u[13] += 1
      else
        u[12] -= 1
      end
      p[7]==4 ? u[49] += 1 : 0
      p[8]==4 ? u[50] += 1 : 0
  elseif temp<cum_rate[21]
      if maximum([p[7],p[8]])>5
        u[13] -= 1
        u[14] += 1
      else
        u[13] -= 1
      end
      p[7]==5 ? u[49] += 1 : 0
      p[8]==5 ? u[50] += 1 : 0
  elseif temp<cum_rate[22]
      if maximum([p[7],p[8]])>6
        u[14] -= 1
        u[15] += 1
      else
        u[14] -= 1
      end
      p[7]==6 ? u[49] += 1 : 0
      p[8]==6 ? u[50] += 1 : 0
  elseif temp<cum_rate[23]
      if maximum([p[7],p[8]])>7
        u[15] -= 1
        u[16] += 1
      else
        u[15] -= 1
      end
      p[7]==7 ? u[49] += 1 : 0
      p[8]==7 ? u[50] += 1 : 0
  elseif temp<cum_rate[24]
      if maximum([p[7],p[8]])>8
        u[16] -= 1
        u[17] += 1
      else
        u[16] -= 1
      end
      p[7]==8 ? u[49] += 1 : 0
      p[8]==8 ? u[50] += 1 : 0
  elseif temp<cum_rate[25]
      if maximum([p[7],p[8]])>9
        u[17] -= 1
        u[18] += 1
      else
        u[17] -= 1
      end
      p[7]==9 ? u[49] += 1 : 0
      p[8]==9 ? u[50] += 1 : 0
  elseif temp<cum_rate[26]
      if maximum([p[7],p[8]])>10
        u[18] -= 1
        u[19] += 1
      else
        u[18] -= 1
      end
      p[7]==10 ? u[49] += 1 : 0
      p[8]==10 ? u[50] += 1 : 0
  elseif temp<cum_rate[27]
      if maximum([p[7],p[8]])>11
        u[19] -= 1
        u[20] += 1
      else
        u[19] -= 1
      end
      p[7]==11 ? u[49] += 1 : 0
      p[8]==11 ? u[50] += 1 : 0
  elseif temp<cum_rate[28]
      if maximum([p[7],p[8]])>12
        u[20] -= 1
        u[21] += 1
      else
        u[20] -= 1
      end
      p[7]==12 ? u[49] += 1 : 0
      p[8]==12 ? u[50] += 1 : 0
  elseif temp<cum_rate[29]
      if maximum([p[7],p[8]])>13
        u[21] -= 1
        u[22] += 1
      else
        u[21] -= 1
      end
      p[7]==13 ? u[49] += 1 : 0
      p[8]==13 ? u[50] += 1 : 0
  elseif temp<cum_rate[30]
      if maximum([p[7],p[8]])>14
        u[22] -= 1
        u[23] += 1
      else
        u[22] -= 1
      end
      p[7]==14 ? u[49] += 1 : 0
      p[8]==14 ? u[50] += 1 : 0
  elseif temp<cum_rate[31]
      if maximum([p[7],p[8]])>15
        u[23] -= 1
        u[24] += 1
      else
        u[23] -= 1
      end
      p[7]==15 ? u[49] += 1 : 0
      p[8]==15 ? u[50] += 1 : 0
  elseif temp<cum_rate[32]
      if maximum([p[7],p[8]])>16
        u[24] -= 1
        u[25] += 1
      else
        u[24] -= 1
      end
      p[7]==16 ? u[49] += 1 : 0
      p[8]==16 ? u[50] += 1 : 0
  elseif temp<cum_rate[33]
      if maximum([p[7],p[8]])>17
        u[25] -= 1
        u[26] += 1
      else
        u[25] -= 1
      end
      p[7]==17 ? u[49] += 1 : 0
      p[8]==17 ? u[50] += 1 : 0
  elseif temp<cum_rate[34]
      if maximum([p[7],p[8]])>18
        u[26] -= 1
        u[27] += 1
      else
        u[26] -= 1
      end
      p[7]==18 ? u[49] += 1 : 0
      p[8]==18 ? u[50] += 1 : 0
  elseif temp<cum_rate[35]
      if maximum([p[7],p[8]])>19
        u[27] -= 1
        u[28] += 1
      else
        u[27] -= 1
      end
      p[7]==19 ? u[49] += 1 : 0
      p[8]==19 ? u[50] += 1 : 0
  elseif temp<cum_rate[36]
      u[28] -= 1
      p[7]==20 ? u[49] += 1 : 0
      p[8]==20 ? u[50] += 1 : 0
  elseif temp<cum_rate[37]
      if maximum([p[7],p[8]])>1
        u[30] += 1
        u[2] -= 1
      else
        u[2] -=1
      end
      p[7]==1 ? u[51] += 1 : 0
      p[8]==1 ? u[52] += 1 : 0
  elseif temp<cum_rate[38]
      if maximum([p[7],p[8]])>2
        u[30] -= 1
        u[31] += 1
      else
        u[30] -= 1
      end
      p[7]==2 ? u[51] += 1 : 0
      p[8]==2 ? u[52] += 1 : 0
  elseif temp<cum_rate[39]
      if maximum([p[7],p[8]])>3
        u[31] -= 1
        u[32] += 1
      else
        u[31] -= 1
      end
      p[7]==3 ? u[51] += 1 : 0
      p[8]==3 ? u[52] += 1 : 0
  elseif temp<cum_rate[40]
      if maximum([p[7],p[8]])>4
        u[32] -= 1
        u[33] += 1
      else
        u[32] -= 1
      end
      p[7]==4 ? u[51] += 1 : 0
      p[8]==4 ? u[52] += 1 : 0
  elseif temp<cum_rate[41]
      if maximum([p[7],p[8]])>5
        u[33] -= 1
        u[34] += 1
      else
        u[33] -= 1
      end
      p[7]==5 ? u[51] += 1 : 0
      p[8]==5 ? u[52] += 1 : 0
  elseif temp<cum_rate[42]
      if maximum([p[7],p[8]])>6
        u[34] -= 1
        u[35] += 1
      else
        u[34] -= 1
      end
      p[7]==6 ? u[51] += 1 : 0
      p[8]==6 ? u[52] += 1 : 0
  elseif temp<cum_rate[43]
      if maximum([p[7],p[8]])>7
        u[35] -= 1
        u[36] += 1
      else
        u[35] -= 1
      end
      p[7]==7 ? u[51] += 1 : 0
      p[8]==7 ? u[52] += 1 : 0
  elseif temp<cum_rate[44]
      if maximum([p[7],p[8]])>8
        u[36] -= 1
        u[37] += 1
      else
        u[36] -= 1
      end
      p[7]==8 ? u[51] += 1 : 0
      p[8]==8 ? u[52] += 1 : 0
  elseif temp<cum_rate[45]
      if maximum([p[7],p[8]])>9
        u[37] -= 1
        u[38] += 1
      else
        u[37] -= 1
      end
      p[7]==9 ? u[51] += 1 : 0
      p[8]==9 ? u[52] += 1 : 0
  elseif temp<cum_rate[46]
      if maximum([p[7],p[8]])>10
        u[38] -= 1
        u[39] += 1
      else
        u[38] -= 1
      end
      p[7]==10 ? u[51] += 1 : 0
      p[8]==10 ? u[52] += 1 : 0
  elseif temp<cum_rate[47]
      if maximum([p[7],p[8]])>11
        u[39] -= 1
        u[40] += 1
      else
        u[39] -= 1
      end
      p[7]==11 ? u[51] += 1 : 0
      p[8]==11 ? u[52] += 1 : 0
  elseif temp<cum_rate[48]
      if maximum([p[7],p[8]])>12
        u[40] -= 1
        u[41] += 1
      else
        u[40] -= 1
      end
      p[7]==12 ? u[51] += 1 : 0
      p[8]==12 ? u[52] += 1 : 0
  elseif temp<cum_rate[49]
      if maximum([p[7],p[8]])>13
        u[41] -= 1
        u[42] += 1
      else
        u[41] -= 1
      end
      p[7]==13 ? u[51] += 1 : 0
      p[8]==13 ? u[52] += 1 : 0
  elseif temp<cum_rate[50]
      if maximum([p[7],p[8]])>14
        u[42] -= 1
        u[43] += 1
      else
        u[42] -= 1
      end
      p[7]==14 ? u[51] += 1 : 0
      p[8]==14 ? u[52] += 1 : 0
  elseif temp<cum_rate[51]
      if maximum([p[7],p[8]])>15
        u[43] -= 1
        u[44] += 1
      else
        u[43] -= 1
      end
      p[7]==15 ? u[51] += 1 : 0
      p[8]==15 ? u[52] += 1 : 0
  elseif temp<cum_rate[52]
      if maximum([p[7],p[8]])>16
        u[44] -= 1
        u[45] += 1
      else
        u[44] -= 1
      end
      p[7]==16 ? u[51] += 1 : 0
      p[8]==16 ? u[52] += 1 : 0
  elseif temp<cum_rate[53]
      if maximum([p[7],p[8]])>17
        u[45] -= 1
        u[46] += 1
      else
        u[45] -= 1
      end
      p[7]==17 ? u[51] += 1 : 0
      p[8]==17 ? u[52] += 1 : 0
  elseif temp<cum_rate[54]
      if maximum([p[7],p[8]])>18
        u[46] -= 1
        u[47] += 1
      else
        u[46] -= 1
      end
      p[7]==18 ? u[51] += 1 : 0
      p[8]==18 ? u[52] += 1 : 0
  elseif temp<cum_rate[55]
      if maximum([p[7],p[8]])>19
        u[47] -= 1
        u[48] += 1
      else
        u[47] -= 1
      end
      p[7]==19 ? u[51] += 1 : 0
      p[8]==19 ? u[52] += 1 : 0
  elseif temp<cum_rate[56]
      u[48] -= 1
      p[7]==20 ? u[51] += 1 : 0
      p[8]==20 ? u[52] += 1 : 0
  end
    t = t+t_step
    push!(out,copy(u))
    push!(t_out,copy(t))
  end
  return out, t_out
end

### Dampening with cXR norepression hypotheis
function sim_delay_damp_cXR_norep(p::Array{Float64,1},tspan::Tuple{Float64,Float64},u0::Array{Int64,1})
  t::Float64=tspan[1]
  i::Int64=0
  out::Array{Array{Int64,1}} = []
  rate_out::Array{Array{Float64,1}} = []
  t_out::Array{Float64,1} = []
  push!(out, u0)
  push!(t_out,t)
  u::Array{Int64,1}=copy(u0)
  #############################################################
  rate::Array{Float64,1}=zeros(56)
  while (t<tspan[2])
    i+=1
    p[24] = 0
    p[29] = 1
    p[25] = 0
    p[30] = 1
    # Xist transcription
    if t<=100
		rate[1] = p[21]*(p[24]+p[29]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11]))
		rate[2] = p[21]*(p[25]+p[30]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11]))
    else
		rate[1] = p[21]*(p[24]+p[29]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[7]^p[13]/(u[7]^p[13]+(p[22]*p[14])^p[13])))
		rate[2] = p[21]*(p[25]+p[30]*(0.5*(u[5]+u[6]))^p[11]/((0.5*(u[5]+u[6]))^p[11]+(p[23]*p[12])^p[11])*(1-u[8]^p[13]/(u[8]^p[13]+(p[22]*p[14])^p[13])))
    end
    # cXR production
    if t<=100
		rate[3] = p[22]
		rate[4] = p[22]
    else
        rate[3] = p[22]*(1-u[51]^p[5]/(u[51]^p[5]+(p[21]*p[6])^p[5]))
        rate[4] = p[22]*(1-u[49]^p[5]/(u[49]^p[5]+(p[21]*p[6])^p[5]))
    end
    # tXA production
    if t<=100
        rate[5] = p[23]
        rate[6] = p[23]
    else
        rate[5] = p[23]*(1-u[52]^p[3]/(u[52]^p[3]+(p[21]*p[4])^p[3]))
        rate[6] = p[23]*(1-u[50]^p[3]/(u[50]^p[3]+(p[21]*p[4])^p[3]))
    end
    # degradation
    rate[7] = 0.1733*u[3]
    rate[8] = 0.1733*u[4]
    rate[9] = 1*u[7]
    rate[10] = 1*u[8]
    rate[11] = 1*u[5]
    rate[12] = 1*u[6]
    rate[13] = 1*u[49]
    rate[14] = 1*u[50]
    rate[15] = 1*u[51]
    rate[16] = 1*u[52]
    ### silencing intermediates on chr1
    if t>100
        rate[17] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[1],0)
        rate[18] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[10],0)
        rate[19] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[11],0)
        rate[20] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[12],0)
        rate[21] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[13],0)
        rate[22] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[14],0)
        rate[23] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[15],0)
        rate[24] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[16],0)
        rate[25] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[17],0)
        rate[26] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[18],0)
        rate[27] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[19],0)
        rate[28] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[20],0)
        rate[29] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[21],0)
        rate[30] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[22],0)
        rate[31] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[23],0)
        rate[32] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[24],0)
        rate[33] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[25],0)
        rate[34] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[26],0)
        rate[35] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[27],0)
        rate[36] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[28],0)
    ### silencing intermediates chr 2
        rate[37] = ifelse(maximum([p[7],p[8]])>=1,p[18]*u[2],0)
        rate[38] = ifelse(maximum([p[7],p[8]])>=2,p[18]*u[30],0)
        rate[39] = ifelse(maximum([p[7],p[8]])>=3,p[18]*u[31],0)
        rate[40] = ifelse(maximum([p[7],p[8]])>=4,p[18]*u[32],0)
        rate[41] = ifelse(maximum([p[7],p[8]])>=5,p[18]*u[33],0)
        rate[42] = ifelse(maximum([p[7],p[8]])>=6,p[18]*u[34],0)
        rate[43] = ifelse(maximum([p[7],p[8]])>=7,p[18]*u[35],0)
        rate[44] = ifelse(maximum([p[7],p[8]])>=8,p[18]*u[36],0)
        rate[45] = ifelse(maximum([p[7],p[8]])>=9,p[18]*u[37],0)
        rate[46] = ifelse(maximum([p[7],p[8]])>=10,p[18]*u[38],0)
        rate[47] = ifelse(maximum([p[7],p[8]])>=11,p[18]*u[39],0)
        rate[48] = ifelse(maximum([p[7],p[8]])>=12,p[18]*u[40],0)
        rate[49] = ifelse(maximum([p[7],p[8]])>=13,p[18]*u[41],0)
        rate[50] = ifelse(maximum([p[7],p[8]])>=14,p[18]*u[42],0)
        rate[51] = ifelse(maximum([p[7],p[8]])>=15,p[18]*u[43],0)
        rate[52] = ifelse(maximum([p[7],p[8]])>=16,p[18]*u[44],0)
        rate[53] = ifelse(maximum([p[7],p[8]])>=17,p[18]*u[45],0)
        rate[54] = ifelse(maximum([p[7],p[8]])>=18,p[18]*u[46],0)
        rate[55] = ifelse(maximum([p[7],p[8]])>=19,p[18]*u[47],0)
        rate[56] = ifelse(maximum([p[7],p[8]])>=20,p[18]*u[48],0)
    end
    
    # Gillespie
    r = rand(2)
    # determine time point of next reaction
    t_step::Float64 = 1/sum(rate)*log(1/(1-r[1]))
    # decide which reaction will occur
    temp::Float64 = r[2]*sum(rate)
    cum_rate = cumsum(rate)
    # Execute reaction
    if temp<cum_rate[1]
      u[1] += 1
      u[3] += 1
        if p[7]==0
        u[49] += 1
      end
      if p[8]==0
        u[50] += 1
      end
  elseif temp<cum_rate[2]
      u[2] += 1
      u[4] += 1
      if p[7]==0
        u[51] += 1
      end
      if p[8]==0
        u[52] += 1
      end
  elseif temp<cum_rate[3]
        u[8] += 1
  elseif temp<cum_rate[4]
      u[7] += 1
  elseif temp<cum_rate[5]
        u[6] += 1
  elseif temp<cum_rate[6]
        u[5] += 1
  elseif temp<cum_rate[7]
      if maximum([p[7],p[8]])==0
        u[1] -= 1
      end
        u[3] -=1
  elseif temp<cum_rate[8]
      if maximum([p[7],p[8]])==0
        u[2] -= 1
      end
        u[4] -=1
  elseif temp<cum_rate[9]
        u[7] -= 1
  elseif temp<cum_rate[10]
        u[8] -= 1
  elseif temp<cum_rate[11]
        u[5] -= 1
  elseif temp<cum_rate[12]
        u[6] -= 1
  elseif temp<cum_rate[13]
        u[49] -= 1
  elseif temp<cum_rate[14]
        u[50] -= 1
  elseif temp<cum_rate[15]
        u[51] -= 1
  elseif temp<cum_rate[16]
        u[52] -= 1
  elseif temp<cum_rate[17]
      if maximum([p[7],p[8]])>1
        u[10] += 1
        u[1] -= 1
      else
        u[1] -= 1
      end
      p[7]==1 ? u[49] += 1 : 0
      p[8]==1 ? u[50] += 1 : 0
  elseif temp<cum_rate[18]
      if maximum([p[7],p[8]])>2
        u[10] -= 1
        u[11] += 1
      else
        u[10] -= 1
      end
      p[7]==2 ? u[49] += 1 : 0
      p[8]==2 ? u[50] += 1 : 0
  elseif temp<cum_rate[19]
      if maximum([p[7],p[8]])>3
        u[11] -= 1
        u[12] += 1
      else
        u[11] -= 1
      end
      p[7]==3 ? u[49] += 1 : 0
      p[8]==3 ? u[50] += 1 : 0
  elseif temp<cum_rate[20]
      if maximum([p[7],p[8]])>4
        u[12] -= 1
        u[13] += 1
      else
        u[12] -= 1
      end
      p[7]==4 ? u[49] += 1 : 0
      p[8]==4 ? u[50] += 1 : 0
  elseif temp<cum_rate[21]
      if maximum([p[7],p[8]])>5
        u[13] -= 1
        u[14] += 1
      else
        u[13] -= 1
      end
      p[7]==5 ? u[49] += 1 : 0
      p[8]==5 ? u[50] += 1 : 0
  elseif temp<cum_rate[22]
      if maximum([p[7],p[8]])>6
        u[14] -= 1
        u[15] += 1
      else
        u[14] -= 1
      end
      p[7]==6 ? u[49] += 1 : 0
      p[8]==6 ? u[50] += 1 : 0
  elseif temp<cum_rate[23]
      if maximum([p[7],p[8]])>7
        u[15] -= 1
        u[16] += 1
      else
        u[15] -= 1
      end
      p[7]==7 ? u[49] += 1 : 0
      p[8]==7 ? u[50] += 1 : 0
  elseif temp<cum_rate[24]
      if maximum([p[7],p[8]])>8
        u[16] -= 1
        u[17] += 1
      else
        u[16] -= 1
      end
      p[7]==8 ? u[49] += 1 : 0
      p[8]==8 ? u[50] += 1 : 0
  elseif temp<cum_rate[25]
      if maximum([p[7],p[8]])>9
        u[17] -= 1
        u[18] += 1
      else
        u[17] -= 1
      end
      p[7]==9 ? u[49] += 1 : 0
      p[8]==9 ? u[50] += 1 : 0
  elseif temp<cum_rate[26]
      if maximum([p[7],p[8]])>10
        u[18] -= 1
        u[19] += 1
      else
        u[18] -= 1
      end
      p[7]==10 ? u[49] += 1 : 0
      p[8]==10 ? u[50] += 1 : 0
  elseif temp<cum_rate[27]
      if maximum([p[7],p[8]])>11
        u[19] -= 1
        u[20] += 1
      else
        u[19] -= 1
      end
      p[7]==11 ? u[49] += 1 : 0
      p[8]==11 ? u[50] += 1 : 0
  elseif temp<cum_rate[28]
      if maximum([p[7],p[8]])>12
        u[20] -= 1
        u[21] += 1
      else
        u[20] -= 1
      end
      p[7]==12 ? u[49] += 1 : 0
      p[8]==12 ? u[50] += 1 : 0
  elseif temp<cum_rate[29]
      if maximum([p[7],p[8]])>13
        u[21] -= 1
        u[22] += 1
      else
        u[21] -= 1
      end
      p[7]==13 ? u[49] += 1 : 0
      p[8]==13 ? u[50] += 1 : 0
  elseif temp<cum_rate[30]
      if maximum([p[7],p[8]])>14
        u[22] -= 1
        u[23] += 1
      else
        u[22] -= 1
      end
      p[7]==14 ? u[49] += 1 : 0
      p[8]==14 ? u[50] += 1 : 0
  elseif temp<cum_rate[31]
      if maximum([p[7],p[8]])>15
        u[23] -= 1
        u[24] += 1
      else
        u[23] -= 1
      end
      p[7]==15 ? u[49] += 1 : 0
      p[8]==15 ? u[50] += 1 : 0
  elseif temp<cum_rate[32]
      if maximum([p[7],p[8]])>16
        u[24] -= 1
        u[25] += 1
      else
        u[24] -= 1
      end
      p[7]==16 ? u[49] += 1 : 0
      p[8]==16 ? u[50] += 1 : 0
  elseif temp<cum_rate[33]
      if maximum([p[7],p[8]])>17
        u[25] -= 1
        u[26] += 1
      else
        u[25] -= 1
      end
      p[7]==17 ? u[49] += 1 : 0
      p[8]==17 ? u[50] += 1 : 0
  elseif temp<cum_rate[34]
      if maximum([p[7],p[8]])>18
        u[26] -= 1
        u[27] += 1
      else
        u[26] -= 1
      end
      p[7]==18 ? u[49] += 1 : 0
      p[8]==18 ? u[50] += 1 : 0
  elseif temp<cum_rate[35]
      if maximum([p[7],p[8]])>19
        u[27] -= 1
        u[28] += 1
      else
        u[27] -= 1
      end
      p[7]==19 ? u[49] += 1 : 0
      p[8]==19 ? u[50] += 1 : 0
  elseif temp<cum_rate[36]
      u[28] -= 1
      p[7]==20 ? u[49] += 1 : 0
      p[8]==20 ? u[50] += 1 : 0
  elseif temp<cum_rate[37]
      if maximum([p[7],p[8]])>1
        u[30] += 1
        u[2] -= 1
      else
        u[2] -=1
      end
      p[7]==1 ? u[51] += 1 : 0
      p[8]==1 ? u[52] += 1 : 0
  elseif temp<cum_rate[38]
      if maximum([p[7],p[8]])>2
        u[30] -= 1
        u[31] += 1
      else
        u[30] -= 1
      end
      p[7]==2 ? u[51] += 1 : 0
      p[8]==2 ? u[52] += 1 : 0
  elseif temp<cum_rate[39]
      if maximum([p[7],p[8]])>3
        u[31] -= 1
        u[32] += 1
      else
        u[31] -= 1
      end
      p[7]==3 ? u[51] += 1 : 0
      p[8]==3 ? u[52] += 1 : 0
  elseif temp<cum_rate[40]
      if maximum([p[7],p[8]])>4
        u[32] -= 1
        u[33] += 1
      else
        u[32] -= 1
      end
      p[7]==4 ? u[51] += 1 : 0
      p[8]==4 ? u[52] += 1 : 0
  elseif temp<cum_rate[41]
      if maximum([p[7],p[8]])>5
        u[33] -= 1
        u[34] += 1
      else
        u[33] -= 1
      end
      p[7]==5 ? u[51] += 1 : 0
      p[8]==5 ? u[52] += 1 : 0
  elseif temp<cum_rate[42]
      if maximum([p[7],p[8]])>6
        u[34] -= 1
        u[35] += 1
      else
        u[34] -= 1
      end
      p[7]==6 ? u[51] += 1 : 0
      p[8]==6 ? u[52] += 1 : 0
  elseif temp<cum_rate[43]
      if maximum([p[7],p[8]])>7
        u[35] -= 1
        u[36] += 1
      else
        u[35] -= 1
      end
      p[7]==7 ? u[51] += 1 : 0
      p[8]==7 ? u[52] += 1 : 0
  elseif temp<cum_rate[44]
      if maximum([p[7],p[8]])>8
        u[36] -= 1
        u[37] += 1
      else
        u[36] -= 1
      end
      p[7]==8 ? u[51] += 1 : 0
      p[8]==8 ? u[52] += 1 : 0
  elseif temp<cum_rate[45]
      if maximum([p[7],p[8]])>9
        u[37] -= 1
        u[38] += 1
      else
        u[37] -= 1
      end
      p[7]==9 ? u[51] += 1 : 0
      p[8]==9 ? u[52] += 1 : 0
  elseif temp<cum_rate[46]
      if maximum([p[7],p[8]])>10
        u[38] -= 1
        u[39] += 1
      else
        u[38] -= 1
      end
      p[7]==10 ? u[51] += 1 : 0
      p[8]==10 ? u[52] += 1 : 0
  elseif temp<cum_rate[47]
      if maximum([p[7],p[8]])>11
        u[39] -= 1
        u[40] += 1
      else
        u[39] -= 1
      end
      p[7]==11 ? u[51] += 1 : 0
      p[8]==11 ? u[52] += 1 : 0
  elseif temp<cum_rate[48]
      if maximum([p[7],p[8]])>12
        u[40] -= 1
        u[41] += 1
      else
        u[40] -= 1
      end
      p[7]==12 ? u[51] += 1 : 0
      p[8]==12 ? u[52] += 1 : 0
  elseif temp<cum_rate[49]
      if maximum([p[7],p[8]])>13
        u[41] -= 1
        u[42] += 1
      else
        u[41] -= 1
      end
      p[7]==13 ? u[51] += 1 : 0
      p[8]==13 ? u[52] += 1 : 0
  elseif temp<cum_rate[50]
      if maximum([p[7],p[8]])>14
        u[42] -= 1
        u[43] += 1
      else
        u[42] -= 1
      end
      p[7]==14 ? u[51] += 1 : 0
      p[8]==14 ? u[52] += 1 : 0
  elseif temp<cum_rate[51]
      if maximum([p[7],p[8]])>15
        u[43] -= 1
        u[44] += 1
      else
        u[43] -= 1
      end
      p[7]==15 ? u[51] += 1 : 0
      p[8]==15 ? u[52] += 1 : 0
  elseif temp<cum_rate[52]
      if maximum([p[7],p[8]])>16
        u[44] -= 1
        u[45] += 1
      else
        u[44] -= 1
      end
      p[7]==16 ? u[51] += 1 : 0
      p[8]==16 ? u[52] += 1 : 0
  elseif temp<cum_rate[53]
      if maximum([p[7],p[8]])>17
        u[45] -= 1
        u[46] += 1
      else
        u[45] -= 1
      end
      p[7]==17 ? u[51] += 1 : 0
      p[8]==17 ? u[52] += 1 : 0
  elseif temp<cum_rate[54]
      if maximum([p[7],p[8]])>18
        u[46] -= 1
        u[47] += 1
      else
        u[46] -= 1
      end
      p[7]==18 ? u[51] += 1 : 0
      p[8]==18 ? u[52] += 1 : 0
  elseif temp<cum_rate[55]
      if maximum([p[7],p[8]])>19
        u[47] -= 1
        u[48] += 1
      else
        u[47] -= 1
      end
      p[7]==19 ? u[51] += 1 : 0
      p[8]==19 ? u[52] += 1 : 0
  elseif temp<cum_rate[56]
      u[48] -= 1
      p[7]==20 ? u[51] += 1 : 0
      p[8]==20 ? u[52] += 1 : 0
  end
    t = t+t_step
    push!(out,copy(u))
    push!(t_out,copy(t))
  end
  return out, t_out
end

# Smoothing the output
using Loess
function smooth_timecourse(t::Array{Float64,1},sol::Array{Int64,2},ti::Array{Float64,1})
    s_smooth::Array{Float64,2}=zeros(Float64,size(sol,1),length(ti))
    s_smooth[:,1] = sol[:,1]
    for n=2:length(ti)
        s_smooth[:,n] = mean(sol[:,(t.<ti[n]).&(t.>=ti[n-1])],2)
    end
    return s_smooth
end

# The following functions take a parameter set and simulate up-regulation (XaXa initiatial cond)
# They give out Xist levels on each of the X chromosomes for each simulated cell

# This function simulates WT cells
function sim_ba_gil(p::Array{Float64,1},nr_cells::Int64, sim_time::Float64)
  ti::Array{Float64,1} = collect(0.0:sim_time)
  out::Array{Float64,2} = Array{Float64}(nr_cells,2*length(ti))
  tspan::Tuple{Float64,Float64} = (0.0,sim_time)
  x0 = zeros(52)
  x0[1:8] = [0,0,0,0,round(p[23]),round(p[23]),round(p[22]),round(p[22])]
  x0_int = convert(Array{Int64},x0)
  for c=1:nr_cells
      temp = sim_delay(p,tspan,x0_int)
      t=temp[2]
      s::Array{Int64,2} = hcat(temp[1]...)
      s_short = smooth_timecourse(t,s,ti)
      rat = log10.((s_short[3,:]+0.1)./(s_short[4,:]+0.1))
      # Train regression model
      model = loess(ti, rat, span=0.2)
      # Apply regression model
      vs = predict(model, ti)
      out[c,1:length(ti)] = s_short[3,:]
      out[c,length(ti)+1:2*length(ti)] = s_short[4,:]
  end
  return out
end

# This function simulates the doxycycline induction experiment
function sim_ba_gil_dox(p::Array{Float64,1},nr_cells::Int64, sim_time::Float64)
  tspan::Tuple{Float64,Float64} = (0.0,sim_time)
  ti::Array{Float64,1} = collect(0.0:sim_time)
  out::Array{Float64,2} = Array{Float64}(nr_cells,2*length(ti))
  x0 = zeros(52)
  x0[1:8] = [0,0,0,0,round(p[23]),round(p[23]),round(p[22]),round(p[22])]
  x0_int = convert(Array{Int64},x0)
  for c=1:nr_cells
      temp = sim_delay_dox(p,tspan,x0_int)
      t=temp[2]
      s::Array{Int64,2} = hcat(temp[1]...)
      s_short = smooth_timecourse(t,s,ti)
      rat = log10.((s_short[3,:]+0.1)./(s_short[4,:]+0.1))
      # Train regression model
      model = loess(ti, rat, span=0.2)
      # Apply regression model
      vs = predict(model, ti)
      out[c,1:length(ti)] = s_short[3,:]
      out[c,length(ti)+1:2*length(ti)] = s_short[4,:]
  end
  return out
end

# This function simulates the artificial biallelic induction experiment
function sim_ba_gil_ba(p::Array{Float64,1},nr_cells::Int64, sim_time::Float64)
  ti::Array{Float64,1} = collect(0.0:sim_time)
  tspan::Tuple{Float64,Float64} = (0.0,sim_time)
  out::Array{Float64,2} = Array{Float64}(nr_cells,4*length(ti))
  x0 = zeros(52)
  x0[1:8] = [0,0,0,0,round(p[23]),round(p[23]),round(p[22]),round(p[22])]
  x0_int = convert(Array{Int64},x0)
  for c=1:nr_cells
      temp = sim_delay_ba(p,tspan,x0_int)
      t1=temp[2]
      p[24] = 0
      p[29] = 1/3
      p[25] = 0
      p[30] = 1
      temp2 = sim_delay_skew(p,tspan,x0_int)
      t2=temp2[2]
      s1::Array{Int64,2} = hcat(temp[1]...)
      s_short1 = smooth_timecourse(t1,s1,ti)
      s2::Array{Int64,2} = hcat(temp2[1]...)
      s_short2 = smooth_timecourse(t2,s2,ti)
      out[c,1:length(ti)] = s_short1[3,:]
      out[c,length(ti)+1:2*length(ti)] = s_short1[4,:]
      out[c,2*length(ti)+1:3*length(ti)] = s_short2[3,:]
      out[c,3*length(ti)+1:4*length(ti)] = s_short2[4,:]
  end
  return out
end

# This function simulates aneuploid cells
function sim_ba_gil_aneu(p::Array{Float64,1},nr_cells::Int64, sim_time::Float64)
  tspan::Tuple{Float64,Float64} = (0.0,sim_time)
  ti::Array{Float64,1} = collect(0.0:sim_time)
  out::Array{Float64,2} = Array{Float64}(nr_cells,4*length(ti))
  if p[26] == 1
      x0 = zeros(52)
      x0[1:8] = [0,0,0,0,round(p[23]),0,round(p[22]),0]
  elseif p[26] == 2
      x0 = zeros(52)
      x0[1:8] = [0,0,0,0,round(p[23]),round(p[23]),round(p[22]),round(p[22])]
  elseif p[26] == 3
      x0 = zeros(102)
      x0[1:8] = [0,0,0,0,round(p[23]),round(p[23]),round(p[22]),round(p[22])]
      x0[53] = 0
      x0[54] = round(p[23])
      x0[55] = round(p[22])
      x0[56:102] = 0
  elseif p[26] == 4
      x0 = zeros(102)
      x0[1:8] = [0,0,0,0,round(p[23]),round(p[23]),round(p[22]),round(p[22])]
      x0[53] = 0
      x0[54] = round(p[23])
      x0[57] = round(p[23])
      x0[55] = round(p[22])
      x0[58] = round(p[22])
      x0[59:102] = 0
  end
  x0_int = convert(Array{Int64},x0)
  for c=1:nr_cells
      temp = sim_delay_aneu(p,tspan,x0_int)
      t=temp[2]
      s::Array{Int64,2} = hcat(temp[1]...)
      s_short = smooth_timecourse(t,s,ti)
      rat = log10.((s_short[3,:]+0.1)./(s_short[4,:]+0.1))
      # Train regression model
      model = loess(ti, rat, span=0.2)
      # Apply regression model
      vs = predict(model, ti)
      out[c,1:length(ti)] = s_short[3,:]
      out[c,length(ti)+1:2*length(ti)] = s_short[4,:]
      if p[26]>=3
          out[c,2*length(ti)+1:3*length(ti)] = s_short[101,:]
          out[c,3*length(ti)+1:4*length(ti)] = s_short[102,:]
      else
          out[c,2*length(ti)+1:3*length(ti)] = 0
          out[c,3*length(ti)+1:4*length(ti)] = 0
      end
  end
  return out
end

# This function simulates polyploid cells with the XA dilution hypothesis
function sim_ba_gil_polypl_dil(p::Array{Float64,1},nr_cells::Int64, sim_time::Float64)
    tspan::Tuple{Float64,Float64} = (0.0,sim_time)
    ti::Array{Float64,1} = collect(0.0:sim_time)
    out::Array{Float64,2} = Array{Float64}(nr_cells,4*length(ti))
  if p[26] == 2
      x0 = zeros(52)
      x0[1:8] = [0,0,0,0,round(p[23]),round(p[23]),round(p[22]),round(p[22])]
  elseif p[26] == 3
      x0 = zeros(102)
      x0[1:8] = [0,0,0,0,round(p[23]),round(p[23]),round(p[22]),round(p[22])]
      x0[53] = 0
      x0[54] = round(p[23])
      x0[55] = round(p[22])
      x0[56:102] = 0
  elseif p[26] == 4
      x0 = zeros(102)
      x0[1:8] = [0,0,0,0,round(p[23]),round(p[23]),round(p[22]),round(p[22])]
      x0[53] = 0
      x0[54] = round(p[23])
      x0[57] = round(p[23])
      x0[55] = round(p[22])
      x0[58] = round(p[22])
      x0[59:102] = 0
  end
  x0_int = convert(Array{Int64},x0)
  for c=1:nr_cells
      temp = sim_delay_polypl_dil(p,tspan,x0_int)
      t=temp[2]
      s::Array{Int64,2} = hcat(temp[1]...)
      s_short = smooth_timecourse(t,s,ti)
      rat = log10.((s_short[3,:]+0.1)./(s_short[4,:]+0.1))
      # Train regression model
      model = loess(ti, rat, span=0.2)
      # Apply regression model
      vs = predict(model, ti)
      out[c,1:length(ti)] = s_short[3,:]
      out[c,length(ti)+1:2*length(ti)] = s_short[4,:]
      if p[26]>=3
          out[c,2*length(ti)+1:3*length(ti)] = s_short[101,:]
          out[c,3*length(ti)+1:4*length(ti)] = s_short[102,:]
      else
          out[c,2*length(ti)+1:3*length(ti)] = 0
          out[c,3*length(ti)+1:4*length(ti)] = 0
      end
  end
  return out
end

# This function simulates polyploid cells with the XA repression hypothesis
function sim_ba_gil_polypl_repXA(p::Array{Float64,1},nr_cells::Int64, sim_time::Float64)
    tspan::Tuple{Float64,Float64} = (0.0,sim_time)
    ti::Array{Float64,1} = collect(0.0:sim_time)
    out::Array{Float64,2} = Array{Float64}(nr_cells,4*length(ti))
  if p[26] == 2
      x0 = zeros(52)
      x0[1:8] = [0,0,0,0,round(p[23]),round(p[23]),round(p[22]),round(p[22])]
  elseif p[26] == 3
      x0 = zeros(102)
      x0[1:8] = [0,0,0,0,round(p[23]),round(p[23]),round(p[22]),round(p[22])]
      x0[53] = 0
      x0[54] = round(p[23])
      x0[55] = round(p[22])
      x0[56:102] = 0
  elseif p[26] == 4
      x0 = zeros(102)
      x0[1:8] = [0,0,0,0,round(p[23]),round(p[23]),round(p[22]),round(p[22])]
      x0[53] = 0
      x0[54] = round(p[23])
      x0[57] = round(p[23])
      x0[55] = round(p[22])
      x0[58] = round(p[22])
      x0[59:102] = 0
  end
  x0_int = convert(Array{Int64},x0)
  for c=1:nr_cells
      temp = sim_delay_polypl_repXA(p,tspan,x0_int)
      t=temp[2]
      s::Array{Int64,2} = hcat(temp[1]...)
      s_short = smooth_timecourse(t,s,ti)
      rat = log10.((s_short[3,:]+0.1)./(s_short[4,:]+0.1))
      # Train regression model
      model = loess(ti, rat, span=0.2)
      # Apply regression model
      vs = predict(model, ti)
      out[c,1:length(ti)] = s_short[3,:]
      out[c,length(ti)+1:2*length(ti)] = s_short[4,:]
      if p[26]>=3
          out[c,2*length(ti)+1:3*length(ti)] = s_short[101,:]
          out[c,3*length(ti)+1:4*length(ti)] = s_short[102,:]
      else
          out[c,2*length(ti)+1:3*length(ti)] = 0
          out[c,3*length(ti)+1:4*length(ti)] = 0
      end
  end
  return out
end

# This function simulates dampening with the cXR upregulation hypothesis
function sim_ba_gil_damp_cXR_upreg(p::Array{Float64,1},nr_cells::Int64, sim_time::Float64)
  ti::Array{Float64,1} = collect(0.0:sim_time)
  tspan::Tuple{Float64,Float64} = (0.0,sim_time)
  out::Array{Float64,2} = Array{Float64}(nr_cells,2*length(ti))
  x0 = zeros(52)
  x0[1:8] = [0,0,0,0,round(p[23]),round(p[23]),0,0]
  x0_int = convert(Array{Int64},x0)
  for c=1:nr_cells
      temp = sim_delay_damp_cXR_upreg(p,tspan,x0_int)
      t=temp[2]
      s::Array{Int64,2} = hcat(temp[1]...)
      s_short = smooth_timecourse(t,s,ti)
      out[c,1:length(ti)] = s_short[3,:]
      out[c,length(ti)+1:2*length(ti)] = s_short[4,:]
  end
  return out
end

# This function simulates dampening with the cXR dampening hypothesis
function sim_ba_gil_damp_cXR_damp(p::Array{Float64,1},nr_cells::Int64, sim_time::Float64)
  ti::Array{Float64,1} = collect(0.0:sim_time)
  tspan::Tuple{Float64,Float64} = (0.0,sim_time)
  out::Array{Float64,2} = Array{Float64}(nr_cells,2*length(ti))
  x0 = zeros(52)
  x0[1:8] = [0,0,0,0,round(p[23]),round(p[23]),round(p[22]),round(p[22])]
  x0_int = convert(Array{Int64},x0)
  for c=1:nr_cells
      temp = sim_delay_damp_cXR_damp(p,tspan,x0_int)
      t=temp[2]
      s::Array{Int64,2} = hcat(temp[1]...)
      s_short = smooth_timecourse(t,s,ti)
      out[c,1:length(ti)] = s_short[3,:]
      out[c,length(ti)+1:2*length(ti)] = s_short[4,:]
  end
  return out
end

# This function simulates dampening with the cXR norepression hypothesis
function sim_ba_gil_damp_cXR_norep(p::Array{Float64,1},nr_cells::Int64, sim_time::Float64)
  ti::Array{Float64,1} = collect(0.0:sim_time)
  tspan::Tuple{Float64,Float64} = (0.0,sim_time)
  out::Array{Float64,2} = Array{Float64}(nr_cells,2*length(ti))
  x0 = zeros(52)
  x0[1:8] = [0,0,0,0,round(p[23]),round(p[23]),round(p[22]),round(p[22])]
  x0_int = convert(Array{Int64},x0)
  for c=1:nr_cells
      temp = sim_delay_damp_cXR_norep(p,tspan,x0_int)
      t=temp[2]
      s::Array{Int64,2} = hcat(temp[1]...)
      s_short = smooth_timecourse(t,s,ti)
      out[c,1:length(ti)] = s_short[3,:]
      out[c,length(ti)+1:2*length(ti)] = s_short[4,:]
  end
  return out
end
