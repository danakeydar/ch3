from  model import *
from conf import *
import dataset as Dataset
from torchmetrics import SpearmanCorrCoef

load_model_ID = 0

train, validation, test, validation_ch3_blind, test_ch3_blind, validation_e_blind, test_e_blind = Dataset.read_data_sets(
    filename_sequence=Conf.filename_sequence,
    filename_expression=Conf.filename_expression,
    filename_labels=Conf.filename_labels,
    filename_dist=Conf.filename_dist,
    train_portion_subjects=Conf.train_portion_probes,
    train_portion_probes=Conf.train_portion_probes, validation_portion_subjects=Conf.validation_portion_subjects,
    validation_portion_probes=Conf.validation_portion_probes, directory='', load_model_ID=load_model_ID)

multi_model = MultiModel()
multi_model = multi_model.float()

lr = 0.001
optimizer = optim.Adam(multi_model.parameters(), lr=lr)
loss_fn = nn.L1Loss()

print(f"train.num_examples= {train.num_examples }")
print(f"Conf.batch_size = {Conf.batch_size}")
print(f"train.num_examples // Conf.batch_size = {train.num_examples // Conf.batch_size}")

train_losses = []; val_losses = []

for epoch in range(Conf.epochs):
  epoch_val_losses=[]; epoch_train_losses=[]; epoch_train_corrs =[]; epoch_val_corrs=[]
  for ii in range(train.num_examples // Conf.batch_size):
    # multi_model.train()
    train_seq_batch, train_exp_batch, train_dist_batch, train_labels_batch = get_next_batch(train, Conf.batch_size)
    optimizer.zero_grad()
    train_predictions = multi_model(train_seq_batch.float(), train_dist_batch.float(), train_exp_batch.float())
    
    batch_train_loss = loss_fn(train_predictions, train_labels_batch)
    epoch_train_losses.append(batch_train_loss.detach().numpy())
    train_spearman = SpearmanCorrCoef()
    batch_train_corr = train_spearman(train_predictions.squeeze(), train_labels_batch.squeeze())
    epoch_train_corrs.append(batch_train_corr)
    
    batch_train_loss.backward()
    optimizer.step()
    
    # multi_model.eval()
    with torch.no_grad():
      val_seq_batch, val_exp_batch, val_dist_batch, val_labels_batch = get_next_batch(validation, Conf.batch_size*2)
      val_predictions = multi_model(val_seq_batch, val_dist_batch, val_exp_batch)
      batch_val_loss = loss_fn(val_predictions, val_labels_batch)
      epoch_val_losses.append(batch_val_loss)
      val_spearman = SpearmanCorrCoef()
      batch_val_corr = val_spearman(val_predictions.squeeze(), val_labels_batch.squeeze())
      epoch_val_corrs.append(batch_val_corr)

  epoch_train_avg_loss = np.mean(epoch_train_losses)
  train_losses.append(epoch_train_avg_loss)

  epoch_val_avg_loss = np.mean(epoch_val_losses)
  val_losses.append(epoch_val_avg_loss)


  print(f"-------------------------------Epoch {epoch +1}/{Conf.epochs}-------------------------------")
  print("Train      | Loss: {:.3f} Corr: {:.3f}".format(epoch_train_avg_loss, np.mean(epoch_train_corrs)))
  print("Validation | Loss: {:.3f} Corr: {:.3f}".format(epoch_val_avg_loss, np.mean(epoch_val_corrs)))
  print("-----------------------------------------------------------------------")

print("Finish training")
print(f"train_losses:{train_losses}")

plt.plot(val_losses, "-b", label="val")
plt.plot(train_losses, "-r", label="train")
plt.legend()
plt.savefig("../res/results_epochs_{epochs}_bs_{batch_size}".format(epochs=Conf.epochs, batch_size=Conf.batch_size))